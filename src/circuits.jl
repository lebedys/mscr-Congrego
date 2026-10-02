#=
Created on 08/12/2022 20:54:46
Last update: 2023-02-06

@author: Michiel Stock
michielfmstock@gmail.com

A basic module to represent the circuits.
=#

export Circuit, PartCatalogue,
    logic_functions, simulate, output2logic,
    n_phenotypes, pheno_index, pheno_label, pheno_labels,
    sequencing_price,
    levenshtein, normalised_levenshtein,
    extract_adjacency!, normalised_graph_distance

import JSON


# ---- Circuit ----

struct Circuit
    parts::Vector{Symbol}
end

# string-based constructor
Circuit(s::AbstractString) = Circuit(Symbol.(split(s, "-")))

Base.length(c::Circuit) = length(c.parts)
Base.get(c::Circuit, i) = c.parts[i]

Base.show(io::IO, c::Circuit) = print(io, "-- " * join(c.parts, " -- ") * " --")


# ---- Part Catalogue ----

# part classification
@enum _PartClass::UInt8 begin
    _EXTERNALLY_ACTIVATABLE_PROMOTER
    _EXTERNALLY_REPRESSIBLE_PROMOTER
    _INTERNALLY_ACTIVATABLE_PROMOTER
    _INTERNALLY_REPRESSIBLE_PROMOTER
    _TRANSCRIPTION_FACTOR
    _TERMINATOR
    _INDUCIBLE_TERMINATOR
    _SPACER
end

struct _PartDescriptor
    class::_PartClass
    is_roadblocker::Bool
    regulator::Symbol     # (:_none if N/A)
end


"""
Collection of parts and their regulatory relationships.
"""
struct PartCatalogue
    parts::Set{Symbol}
    inputs::Vector{Symbol}
    activatable_promoters::Dict{Symbol,Symbol} # :promoter => :activatingTF
    repressible_promoters::Dict{Symbol,Symbol} # :promoter => :repressingTF
    terminators::Set{Symbol}
    roadblockers::Set{Symbol}
    inducible_terminators::Dict{Symbol,Symbol} # :terminator => :TF
    part_length_map::Dict{Symbol,Int} # :part => length

    # derived fields
    promoters::Set{Symbol}
    transcription_factors::Set{Symbol}
    functional_parts::Set{Symbol}
    part_descriptors::Dict{Symbol, _PartDescriptor}
    node_order::Vector{Symbol}
    node_index::Dict{Symbol, Int}

    function PartCatalogue(; parts::Set{Symbol},
        inputs::Vector{Symbol},
        activatable_promoters::Dict{Symbol,Symbol},
        repressible_promoters::Dict{Symbol,Symbol},
        terminators::Set{Symbol},
        roadblockers::Set{Symbol},
        inducible_terminators::Dict{Symbol,Symbol},
        part_length_map::Dict{Symbol,Int})

        @assert isempty(parts ∩ Set(inputs)) "parts and inputs must not overlap"
        @assert keys(activatable_promoters) ⊆ parts "activatable promoter keys must be in parts"
        @assert values(activatable_promoters) ⊆ (parts ∪ Set(inputs)) "activatable promoter transcription factors must be in parts or inputs"
        @assert keys(repressible_promoters) ⊆ parts "repressible promoter keys must be in parts"
        @assert values(repressible_promoters) ⊆ (parts ∪ Set(inputs)) "repressible promoter transcription factors must be in parts or inputs"
        @assert terminators ⊆ parts "terminators must be in parts"
        @assert roadblockers ⊆ parts "roadblockers must be in parts"
        @assert keys(inducible_terminators) ⊆ parts "inducible terminators must be in parts"
        @assert values(inducible_terminators) ⊆ parts "inducible terminator transcription factors must be in parts or inputs"
        @assert parts ⊆ keys(part_length_map) "all parts must have lengths in part_length_map"

        promoters = Set{Symbol}( keys(activatable_promoters) ∪ keys(repressible_promoters) )

        # set of unique TFs that induce promoters/terminators
        all_values = union(
            values(activatable_promoters),
            values(repressible_promoters),
            values(inducible_terminators),
        )
        # exclude inputs (only include TFs available in parts)
        transcription_factors = Set{Symbol}(filter(v -> v ∈ parts, all_values))

        # precompute part descriptors to avoid the runtime hit of dictionary lookups during simulate()
        input_syms_set = Set{Symbol}(inputs)
        part_descriptors = Dict{Symbol, _PartDescriptor}()
        for part in parts
            is_roadblocker = part ∈ roadblockers
            if part ∈ keys(activatable_promoters)
                activator = activatable_promoters[part]
                part_class = activator ∈ input_syms_set ? _EXTERNALLY_ACTIVATABLE_PROMOTER : _INTERNALLY_ACTIVATABLE_PROMOTER
                part_descriptors[part] = _PartDescriptor(part_class, is_roadblocker, activator)
            elseif part ∈ keys(repressible_promoters)
                repressor = repressible_promoters[part]
                part_class = repressor ∈ input_syms_set ? _EXTERNALLY_REPRESSIBLE_PROMOTER : _INTERNALLY_REPRESSIBLE_PROMOTER
                part_descriptors[part] = _PartDescriptor(part_class, is_roadblocker, repressor)
            elseif part ∈ transcription_factors
                part_descriptors[part] = _PartDescriptor(_TRANSCRIPTION_FACTOR, false, :_none)
            elseif part ∈ terminators
                part_descriptors[part] = _PartDescriptor(_TERMINATOR, false, :_none)
            elseif part ∈ keys(inducible_terminators)
                part_descriptors[part] = _PartDescriptor(_INDUCIBLE_TERMINATOR, false, inducible_terminators[part])
            else
                part_descriptors[part] = _PartDescriptor(_SPACER, false, :_none)
            end
        end

        # all parts that affect circuit behavior (everything except spacers)
        functional_parts = transcription_factors ∪ keys(activatable_promoters) ∪
                           keys(repressible_promoters) ∪ terminators ∪
                           keys(inducible_terminators)

        # node ordering for adjacency matrices: sorted (inputs ∪ TFs), then :output
        node_order = sort(collect(transcription_factors ∪ Set(inputs)))
        push!(node_order, :output)
        node_index = Dict(sym => i for (i, sym) in enumerate(node_order))

        new(parts, inputs,
            activatable_promoters, repressible_promoters,
            terminators, roadblockers, inducible_terminators,
            part_length_map,
            promoters, transcription_factors, functional_parts,
            part_descriptors, node_order, node_index)
    end
end


const single_input_logic_functions = [
    Symbol(false), # FF
    :identity,     # FT
    :not,          # TF
    Symbol(true)   # TT
]

# TODO - rename to include two_input_ in the name
const logic_functions = [
    Symbol(false), # FFFF
    :and,          # FFFT
    :animplyb,     # FFTF
    :a,            # FFTT
    :bnimplya,     # FTFF
    :b,            # FTFT
    :xor,          # FTTF
    :or,           # FTTT
    :nor,          # TFFF
    :xnor,         # TFFT          
    :notb,         # TFTF 
    :bimplya,      # TFTT    
    :nota,         # TTFF 
    :aimplyb,      # TTFT    
    :nand,         # TTTF 
    Symbol(true)   # TTTT       
]
output2logic((s1,s2,s3,s4)) = logic_functions[8s1+4s2+2s3+s4+1]

"For k inputs, return the number of phenotypes (i.e. logical functions) that can be encoded."
n_phenotypes(k) = 2^(2^k)


"""
Return an integer index corresponding to the binary value of the outputs of a given phenotype.

outputs = (true, true, false, false)
        => 1100 => 1*8 + 1*4 + 0*2 + 0*1 + 1
        => 12 + 1 => 13

This is so we can index an array of phenotype labels, e.g. logic_functions above
"""
pheno_index(outputs) = sum(o * 2^(length(outputs)-i) for (i, o) in enumerate(outputs)) + 1

"Return a symbolic label for a phenotype."
function pheno_label(outputs)
    k = Int(log2(length(outputs)))
    # For k=1,2 inputs, we have names for the logic functions they can encode
    k == 1 && return single_input_logic_functions[pheno_index(outputs)]
    k == 2 && return logic_functions[pheno_index(outputs)]
    # Above 2 inputs, not all functions have a name.
    # The label is just the binary representation of the function's outputs under all input conditions.
    # TODO - can we return a logical function of the inputs in Boolean algebraic notation? e.g. (A+B)C
    return Symbol(join(Int.(outputs)))
end

"""
Return a vector of all phenotype labels for k inputs
"""
function pheno_labels(k::Integer)
    # For k=1,2 inputs, we have names for the logic functions they can encode
    k == 1 && return single_input_logic_functions
    k == 2 && return logic_functions
    # Above 2 inputs, not all functions have a name.
    # The label is just the binary representation of the function's outputs under all input conditions.

    return vec([pheno_label(reverse(c)) for c in Iterators.product(fill([false, true], 2^k)...)])
end

"""
Price in GBP for synthesising a part with length `l`

from https://eu.idtdna.com/pages/products/genes-and-gene-fragments/double-stranded-dna-fragments/gblocks-gene-fragments
"""
function sequencing_price(l)
    length_price = [
        (250,  37.0), 
        (500,  37.0),
        (750,  54.0), 
        (1000, 71.0),  
        (1250, 88.0), 
        (1500, 105.0),
        (1750, 125.0), 
        (2000, 150.0), 
        (2250, 175.0),
        (2500, 200.0), 
        (2750, 225.0), 
        (3000, 245.0)]
    i = findfirst(≥((l,0)), length_price)
    isnothing(i) && return Inf
    return length_price[i][2]
end


# ---- Simulate ----

"""
    simulate(circuit::Circuit, inputs::Pair{Symbol,Bool}...; pc::PartCatalogue, verbose=false)

Simulate a circuit with named inputs.
External inputs (e.g. IPTG, aTC) are specified as `Symbol => Bool` pairs.
Unspecified inputs default to false.

# Examples
```julia
simulate(Circuit([:pTac, :tfSp]), :IPTG => true, :aTC => false; pc=pc_m)
simulate(Circuit([:pTac, :tfSp]), :IPTG => true; pc=pc_m)
```
"""
function simulate(circuit::Circuit, inputs::Pair{Symbol,Bool}...; pc::PartCatalogue, verbose=false)
    n = length(circuit)
    state = zeros(Bool, n)
    visited_states = typeof(state)[]

    # get part descriptors for this circuit
    parts = [pc.part_descriptors[part] for part in circuit.parts]

    verbose && @show inputs

    # induce input-driven promoters
    for (i, part) in enumerate(parts)
        if part.class == _EXTERNALLY_ACTIVATABLE_PROMOTER
            # activate if input is true
            state[i] = any(input_val && input_symbol == part.regulator for (input_symbol, input_val) in inputs)
        elseif part.class == _EXTERNALLY_REPRESSIBLE_PROMOTER
            # repress if input is true
            state[i] = !any(input_val && input_symbol == part.regulator for (input_symbol, input_val) in inputs)
        end
    end

    gfp = Bool[]
    while true
        transcr = false
        for (i, part) in enumerate(parts)
            if part.class == _EXTERNALLY_ACTIVATABLE_PROMOTER
                # state set during initialization
                if part.is_roadblocker
                    transcr = state[i]
                else
                    transcr |= state[i]
                end
                verbose && @show part.class, circuit.parts[i], state[i], transcr

            elseif part.class == _EXTERNALLY_REPRESSIBLE_PROMOTER
                # state set during initialization
                if part.is_roadblocker
                    transcr = state[i]
                else
                    transcr |= state[i]
                end
                verbose && @show part.class, circuit.parts[i], state[i], transcr

            elseif part.class == _INTERNALLY_ACTIVATABLE_PROMOTER
                # activate if induced
                state[i] = any(s && p == part.regulator for (s, p) in zip(state, circuit.parts))
                if part.is_roadblocker
                    transcr = state[i]
                else
                    transcr |= state[i]
                end
                verbose && @show part.class, circuit.parts[i], state[i], transcr

            elseif part.class == _INTERNALLY_REPRESSIBLE_PROMOTER
                # repress if induced
                state[i] = !any(s && p == part.regulator for (s, p) in zip(state, circuit.parts))
                if part.is_roadblocker
                    transcr = state[i]
                else
                    transcr |= state[i]
                end
                verbose && @show part.class, circuit.parts[i], state[i], transcr

            elseif part.class == _TRANSCRIPTION_FACTOR
                # transcribe if upstream promoter is active
                state[i] = transcr
                verbose && @show part.class, circuit.parts[i], state[i]

            elseif part.class == _SPACER
                # spacer does nothing
                verbose && @show part.class, circuit.parts[i]

            elseif part.class == _INDUCIBLE_TERMINATOR
                # halt transcription if induced
                if transcr
                    transcr &= !any(s && p == part.regulator for (s, p) in zip(state, circuit.parts))
                end
                verbose && @show part.class, circuit.parts[i], transcr

            elseif part.class == _TERMINATOR
                # halt transcription
                transcr = false
                verbose && @show part.class, circuit.parts[i], transcr
            end
        end

        # check convergence
        if state ∈ visited_states
            if state == last(visited_states)
                # if state is stable, return last activity
                return (transcr, copy(state))
            else
                ind = findfirst(isequal(state), visited_states)
                activation = @view gfp[ind:end]
                all(activation) && return (true, copy(state))
                !any(activation) && return (false, copy(state))
                return (missing, copy(state))
            end
        else
            push!(visited_states, copy(state))
            push!(gfp, transcr)
        end
        verbose && @show state, gfp
    end
end


# Backward-compatible 2-input convenience method
simulate(circuit::Circuit, IPTG::Bool, aTC::Bool; kw...) =
    simulate(circuit, :IPTG => IPTG, :aTC => aTC; kw...)


"""
    simulate(circuit::Circuit, input_syms::AbstractVector{Symbol}; pc::PartCatalogue)

Simulate a circuit for all combinations of the given input symbols.
Returns an NTuple of `2^k` results where k is the number of inputs.

The order of input vector will determine the order of input permutations.
For example, the vector:
    input_syms = [:IPTG,:aTC]

Will yield the truth table:
    | :IPTG | :aTC  |
    | false | false |
    | false | true  |
    | true  | false |
    | true  | true  |

Therefore the returned phenotype will be:
    [output(IPTG=false,aTC=false), output(IPTG=false,aTC=true), output(IPTG=true,aTC=false), output(IPTG=true,aTC=true)]

# Examples
```julia
simulate(Circuit([:pTac, :tfSp]), [:IPTG, :aTC]; pc=pc_m)
simulate(Circuit([:pTac, :tfSp]), [:IPTG]; pc=pc_m)
```
"""
function simulate(circuit::Circuit, input_syms::AbstractVector{Symbol}; pc::PartCatalogue)
    k = length(input_syms)
    combos = Iterators.product(fill([false, true], k)...)
    # reverse(c) so that last input varies fastest (standard truth table order),
    # since Iterators.product varies the first dimension fastest (column-major)
    results = [simulate(circuit, (input_syms .=> reverse(c))...; pc) for c in combos] |> vec
    outputs = NTuple{2^k,Union{Missing,Bool}}(first(r) for r in results)
    states = NTuple{2^k,Vector{Bool}}(last(r) for r in results)
    return (outputs, states)
end


"""
    simulate(circuit::AbstractString, input_syms::AbstractVector{Symbol}; pc::PartCatalogue)

Simulate a circuit represented as a string with parts separated by `-`.

# Examples
```julia
simulate("pTac-tfSp", [:IPTG, :aTC]; pc=pc_m)
```
"""
simulate(circuit::AbstractString, input_syms::AbstractVector{Symbol}; pc::PartCatalogue) =
    Symbol.(split(circuit, "-")) |> Circuit |> c -> simulate(c, input_syms; pc)


# Original 2-input simulate archived for testing against the refactored version
function _simulate_v1(circuit::Circuit, IPTG, aTC; pc::PartCatalogue, verbose=false)
    n = length(circuit)
    state = zeros(Bool, n)
    visited_states = typeof(state)[]
    # add activators
    IPTG && @.(state[:pTac==circuit.parts] = true)
    aTC && @.(state[:pTet==circuit.parts] = true)
    verbose && @show IPTG aTC
    gfp = Bool[]
    while true
        # update the states
        transcr = false
        for (i, part) in enumerate(circuit.parts)
            # roadblocking stops transcription
            if part == :pTac || part == :pTet
                transcr |= state[i]
            elseif part ∈ keys(pc.activatable_promoters)
                # check if there are any active activators
                activator = pc.activatable_promoters[part]
                state[i] = any(s && p == activator for (s, p) in zip(state, circuit.parts))
                transcr |= state[i]
            elseif part ∈ keys(pc.repressible_promoters)
                repressor = pc.repressible_promoters[part]
                state[i] = !any(s && p == repressor for (s, p) in zip(state, circuit.parts))
                if part ∉ pc.roadblockers
                    transcr |= state[i]
                else
                    transcr = state[i]
                end
            elseif part ∈ pc.transcription_factors
                state[i] = transcr
            elseif part ∈ pc.terminators
                transcr = false
            elseif transcr && part ∈ keys(pc.inducible_terminators)
                inducer = pc.inducible_terminators[part]
                transcr &= !any(s && p == inducer for (s, p) in zip(state, circuit.parts))
            end
        end
        if state ∈ visited_states
            # if state is stable, return last activity
            if state == last(visited_states)
                return transcr
            else
                ind = findfirst(isequal(state), visited_states)
                activation = @view gfp[ind:end]
                all(activation) && return true
                !any(activation) && return false
                return missing
            end
        else
            push!(visited_states, copy(state))
            push!(gfp, transcr)
        end
        verbose && @show state, gfp
    end
end


# ----- Edit Distance ----

"""
    levenshtein(a, b)

DP implementation of Levenshten edit distance for two circuits.
- https://www.baeldung.com/cs/levenshtein-distance-computation
"""
function levenshtein(a, b; verbose=false)
    m, n = length(a), length(b)

    # choose smaller circuit for matrix columns
    if m < n; a, b = b, a; m, n = n, m; end
    verbose && @info "levenshtein" a b m n

    # instead of the whole (m,n) matrix,
    # we only keep two rows in memory:
    # the current and the previous
    prev = collect(0:n)
    curr = similar(prev)

    # verbose && @info "row 0 (base)" prev
    verbose && (rows = [copy(prev)])

    for i in 1:m
        curr[1] = i
        for j in 1:n
            cost = a[i] == b[j] ? 0 : 1
            curr[j+1] = min(
                prev[j+1] + 1,
                curr[j] + 1,
                prev[j] + cost
            )
        end
        # verbose && @info "row $i" a_i=a[i] curr=copy(curr)
        verbose && push!(rows, copy(curr))

        # swap row references
        prev, curr = curr, prev
    end

    if verbose
        _print_levenshtein_matrix(a, b, rows)
        @info "result" distance=prev[n+1]
    end

    # after final prev/curr reference swap, 
    # prev refers to the final row in the matrix
    # take its last element as the distance metric
    return prev[n+1]
end

"""
    normalised_levenshtein(a, b)
Edit distance normalised to [0,1]
"""
function normalised_levenshtein(a, b; verbose=false)
    maxlen = max(length(a), length(b))
    maxlen == 0 && return 0.0
    distance = levenshtein(a, b; verbose)
    normed = distance / maxlen

    verbose && @info "normalised_levenshtein" distance maxlen normed
    return normed
end


"""
    extract_adjacency!(adj, parts; pc) → adj

Fill a pre-allocated N×N `Matrix{Int8}` with the signed regulatory adjacency matrix for a circuit. 

Entries ∈ {-1, 0, +1}: activation (+1), repression (-1), none (0). 
This assumes that each regulator either activates or represses, never both.
"""
function extract_adjacency!(adj::Matrix{Int8}, parts::Vector{Symbol}; pc::PartCatalogue, verbose::Bool=false)

    # pre-computed indeces in the adjacency matrix for each node
    node_index = pc.node_index

    # initialise empty
    fill!(adj, Int8(0))

    # Keep track of regulators driving upstream promoters
    # e.g. for the circuit pSrpR-pTet-AraC
    #      the upstream_regulators for AraC: [(SrpR, -1), (aTc, +1)]
    # 
    # Terminators reset this list
    # e.g. for the circuit pTet-T10-AraC
    #      the upstream_regulators for AraC: []
    # 
    # Roadblockers also reset this list
    # e.g. for the circuit pTet-pSrpR-AraC
    #      the upstream_regulators for AraC: [(SrpR, -1)]
    # Because AraC transcription now depends on: 
    #   aTc AND NOT(SrpR)  OR  NOT(aTc) AND NOT(SrpR) 
    # Which logically simplifies to:
    #   NOT(SrpR)
    upstream_regulators = Tuple{Int,Int8}[]  # (regulator_node, regulation_sign)

    for part in parts
        part_descriptor = pc.part_descriptors[part]

        if part_descriptor.class ∈ (_EXTERNALLY_ACTIVATABLE_PROMOTER, _INTERNALLY_ACTIVATABLE_PROMOTER)
            regulator_node = node_index[part_descriptor.regulator]
            part_descriptor.is_roadblocker && empty!(upstream_regulators)
            push!(upstream_regulators, (regulator_node, Int8(1)))

        elseif part_descriptor.class ∈ (_EXTERNALLY_REPRESSIBLE_PROMOTER, _INTERNALLY_REPRESSIBLE_PROMOTER)
            regulator_node = node_index[part_descriptor.regulator]
            part_descriptor.is_roadblocker && empty!(upstream_regulators)
            push!(upstream_regulators, (regulator_node, Int8(-1)))

        elseif part_descriptor.class == _TRANSCRIPTION_FACTOR
            target_node = node_index[part]
            @inbounds for (regulator_node, regulation_sign) in upstream_regulators
                adj[regulator_node, target_node] = regulation_sign
            end

        elseif part_descriptor.class == _TERMINATOR
            empty!(upstream_regulators)

        elseif part_descriptor.class == _INDUCIBLE_TERMINATOR
            if !isempty(upstream_regulators)
                regulator_node = node_index[part_descriptor.regulator]
                push!(upstream_regulators, (regulator_node, Int8(-1)))
            end

        # _SPACER: no effect
        end
    end

    # Output edges
    target_node = node_index[:output]
    @inbounds for (regulator_node, regulation_sign) in upstream_regulators
        adj[regulator_node, target_node] = regulation_sign
    end

    verbose && _print_adj_matrix(adj, pc)
    return adj
end



"""
    normalised_graph_distance(adjacency_a, adjacency_b) → Float64

Hamming distance over active positions of two signed adjacency matrices.
Active positions are those where at least one matrix has a non-zero entry.

Normalised to [0, 1].
"""
function normalised_graph_distance(adjacency_a::Matrix{Int8}, adjacency_b::Matrix{Int8})
    n_active_positions = 0
    n_mismatches = 0
    for i in eachindex(adjacency_a, adjacency_b)
        edge_a, edge_b = adjacency_a[i], adjacency_b[i]
        (edge_a == 0 && edge_b == 0) && continue
        n_active_positions += 1
        (edge_a != edge_b) && (n_mismatches += 1)
    end
    n_active_positions == 0 && return 0.0
    return n_mismatches / n_active_positions
end



# ----- Serialisation & Debugging Utilities ----


# JSON serialisation to allow JSON.print(pc)
function JSON.lower(pc::PartCatalogue)
    _sort(s::AbstractSet)  = sort(collect(s), by=string)
    _sort(d::AbstractDict) = Dict(k => d[k] for k in sort(collect(keys(d)), by=string))
    Dict(
        :parts                 => _sort(pc.parts),
        :inputs                => pc.inputs,
        :activatable_promoters => _sort(pc.activatable_promoters),
        :repressible_promoters => _sort(pc.repressible_promoters),
        :terminators           => _sort(pc.terminators),
        :roadblockers          => _sort(pc.roadblockers),
        :inducible_terminators => _sort(pc.inducible_terminators),
        :promoters             => _sort(pc.promoters),
        :transcription_factors => _sort(pc.transcription_factors),
        :part_length_map       => _sort(pc.part_length_map),
    )
end

# Deserialisation: reconstruct a PartCatalogue from JSON.parsefile output
function PartCatalogue(d::AbstractDict)
    PartCatalogue(
        parts                 = Set(Symbol.(d["parts"])),
        inputs                = Symbol.(d["inputs"]),
        activatable_promoters = Dict(Symbol(k) => Symbol(v) for (k,v) in d["activatable_promoters"]),
        repressible_promoters = Dict(Symbol(k) => Symbol(v) for (k,v) in d["repressible_promoters"]),
        terminators           = Set(Symbol.(d["terminators"])),
        roadblockers          = Set(Symbol.(d["roadblockers"])),
        inducible_terminators = Dict(Symbol(k) => Symbol(v) for (k,v) in d["inducible_terminators"]),
        part_length_map       = Dict(Symbol(k) => Int(v) for (k,v) in d["part_length_map"]),
    )
end

"""
Print DP matrix for debugging
"""
function _print_levenshtein_matrix(a, b, rows)
    buf = IOBuffer()
    syms = vcat(a, b)
    w = isempty(syms) ? 4 : max(4, maximum(length ∘ string, syms) + 1)
    print(buf, lpad("", w), lpad("_", w))
    for bj in b; print(buf, lpad(string(bj), w)); end
    println(buf)
    for (idx, row) in enumerate(rows)
        label = idx == 1 ? "_" : string(a[idx-1])
        print(buf, lpad(label, w))
        for val in row; print(buf, lpad(string(val), w)); end
        println(buf)
    end
    @info "DP matrix\n" * String(take!(buf))
end


"""
Print the signed adjacency matrix for debugging (rows = source, cols = target).
"""
function _print_adj_matrix(adj::Matrix{Int8}, pc::PartCatalogue)
    nodes = pc.node_order
    buf = IOBuffer()
    w = max(7, maximum(length ∘ string, nodes) + 1)
    print(buf, lpad("", w))                          # corner
    for n in nodes; print(buf, lpad(string(n), w)); end
    println(buf)
    for (i, src) in enumerate(nodes)                 # one row per source node
        print(buf, lpad(string(src), w))
        for j in eachindex(nodes)
            v = adj[i, j]
            print(buf, lpad(v == 0 ? "." : string(Int(v)), w))   # 0 → "." so edges pop
        end
        println(buf)
    end
    @info "Adjacency matrix (rows=source, cols=target)\n" * String(take!(buf))
end