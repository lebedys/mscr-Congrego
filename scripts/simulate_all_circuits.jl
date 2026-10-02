#=
Created on 10/01/2023 14:58:12
Last update: 2023-01-30

@author: Michiel Stock
michielfmstock@gmail.com

Simulate all the basic circuits with the parts.
=#

using DrWatson
@quickactivate "Congrego"
using Congrego, Combinatorics, DataFrames, CSV, SparseArrays

using Revise
Revise.revise()


for n_tiles in 1:6

    println("Simulating all $(length(parts)^n_tiles) circuits for $n_tiles tiles...")

    generator = Iterators.product([parts for _ in 1:n_tiles]...)

    n_circuits = length(generator)

    @time results = [(ps, first(simulate(Circuit([ps...]), external_inputs)))
        for ps in generator] |> vec

    circuits = first.(results)
    outputs = last.(results)

    gp_map = Dict(c=>o for (c, o) in zip(circuits, outputs)
                    if !any(ismissing,o))

    #BSON.@save datadir("gp_map.bson") gp_map circuits outputs

    # some basic statistics about the functions

    unique_functions = filter(o->!any(ismissing, o), unique(outputs))

    n_circuits = length(circuits)
    n_circuits_oscillator = count(o->any(ismissing, o), outputs)
    n_circuits_nill = count(o->all(!ismissing,o) && !any(o), outputs)

    circuit_type_counts = Dict(f=>count(isequal(f), outputs) for f in unique_functions)

    df = DataFrame([(tiles..., outp...,
                    any(ismissing, outp))
        for (tiles, outp) in results])

    rename!(df, [("t$i" for i in 1:n_tiles)..., "none", "aTC", "IPTG", "both", "cyclic"])

    CSV.write(datadir("circuits_$(n_tiles)_phenotypes.csv"), df)
end

# can we make a genotype network and assess its stucture?

function make_genotype_network(gp_map, parts)
    n = length(gp_map)
    N = spzeros(UInt, n, n)
    checked = Set{eltype(keys(gp_map))}()
    c = 0
    for circuit in keys(gp_map)
        c += 1
        circuit ∈ checked && continue
        push!(checked, circuit)
        ind, phenotype = gp_map[circuit]
        # only consider non-osc circuits
        for i in 1:length(circuit)
            for part in parts
                neighbor = NTuple{6,Symbol}(i==j ? part : p for (j,p) in enumerate(circuit))
                neighbor ∈ checked && continue
                ind_n, pheno_n = gp_map[neighbor]
                if pheno_n == phenotype
                     N[ind, ind_n] = phenotype
                     N[ind_n, ind] = phenotype
                end
            end
        end
        @show c, circuit
        c ≥ 100 && break
    end
    return N
end

# note: this works extremely slow (but I don't see why)?

#gp_map = Dict(g=>(i,hash(p)) for (i, (g, p)) in enumerate(zip(circuits, outputs)))

#genotype_network = make_genotype_network(gp_map, parts)
