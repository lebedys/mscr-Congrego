
export Generator, CGenerator, number_of_circuits, number_of_parts,
        number_of_slots, gen2parts, part_length_bp, get_price, mean_part_bp,
        MAX_PART_BP, iter, get_parts,
        print_phenotypes, characterize, characterize_with_diversity, summarize,
        phenotype_distribution, plot_phenotypes


using UnicodePlots: barplot
using Plots: bar, plot!, annotate!, text, mm, @colorant_str
# Generator structure
# -------------------

struct Generator{N,T}
    slots::NTuple{N,Set{T}}
end

# define parts in generator explicitly, as a vector
Generator(g::Vector) = Generator(NTuple{length(g),Set{eltype(g[1])}}(Set.(g)))
# random generator of singleton parts
Generator(n, parts) = Generator(NTuple{n,Set{eltype(parts)}}([Set([rand(parts)]) for i in 1:n]))
# random generator of composite parts
CGenerator(n, parts) = Generator(NTuple{n,Set{Vector{eltype(parts)}}}([Set([[rand(parts)]]) for i in 1:n]))


# part types to dispatch on Generator vs CGenerator
abstract type PartType end
struct SingletonParts <: PartType end
struct CompositeParts <: PartType end

part_type(::Generator) = SingletonParts
part_type(::Generator{<:Any,<:AbstractVector}) = CompositeParts

# versions with singleton parts
Base.rand(::Type{<:SingletonParts}, g::Generator) = [rand(slot) for slot in g.slots if !isempty(slot)]
iter(::Type{<:SingletonParts}, g::Generator) = Iterators.product(filter(!isempty, g.slots)...)
Base.collect(::Type{<:SingletonParts}, g::Generator) = iter(g) |> collect |> vec

# versions with composite parts
Base.rand(::Type{<:CompositeParts}, g::Generator) = vcat([rand(slot) for slot in g.slots]...)
# flatten composite slot tuples into a single Vector{Symbol} per circuit
iter(::Type{<:CompositeParts}, g::Generator) = (collect(Symbol, Iterators.flatten(c)) for c in Iterators.product(filter(!isempty, g.slots)...))
Base.collect(::Type{<:CompositeParts}, g::Generator) = iter(g) |> collect |> vec

Base.rand(g::Generator) = rand(part_type(g), g)
iter(g::Generator) = iter(part_type(g), g)
Base.collect(g::Generator) = collect(part_type(g), g)

number_of_circuits(g::Generator) = prod(s->max(1,length(s)), g.slots)
number_of_slots(g::Generator) = count(!isempty, g.slots)
number_of_parts(g::Generator) = sum(length, g.slots)


# TODO adomasl - remove this: no longer necessary now that we use symbols instead of indices for generators
gen2parts(g::Generator) = gen2parts(part_type(g), g)
gen2parts(::Type{<:SingletonParts}, g::Generator) = [collect(slot) for slot in g.slots] |> vec
gen2parts(::Type{<:CompositeParts}, g::Generator) = [collect(slot) for slot in g.slots] |> vec

get_parts(g::Generator) = get_parts(part_type(g), g)
get_parts(::Type{<:SingletonParts}, g::Generator) = collect(reduce(union, g.slots)) |> sort!
get_parts(::Type{<:CompositeParts}, g::Generator) = reduce(union, reduce(union, g.slots)) |> sort!
Base.isempty(g::Generator) = all(isempty, g.slots)

# Base pair length of a part (singleton or composite)
part_length_bp(s::Symbol, pc::PartCatalogue) = pc.part_length_map[s]
part_length_bp(c::AbstractVector{Symbol}, pc::PartCatalogue) = sum(s -> pc.part_length_map[s], c)

# Synthesis price (composes with part_length_bp)
get_price(s::Symbol, pc::PartCatalogue) = sequencing_price(part_length_bp(s, pc))
get_price(c::AbstractVector{Symbol}, pc::PartCatalogue) = sequencing_price(part_length_bp(c, pc))
# of a generator: same parts at different positions!
get_price(g::Generator, pc::PartCatalogue) = sum(mapreduce(s->get_price(s, pc), +, slot, init=0.0) for slot in g.slots)

const MAX_PART_BP = 3000

"""
    mean_part_bp(g::Generator, pc::PartCatalogue)

Mean base pair length across all parts in the generator.
"""
function mean_part_bp(g::Generator, pc::PartCatalogue)
    total = 0
    n = 0
    for slot in g.slots, part in slot
        total += part_length_bp(part, pc)
        n += 1
    end
    return n > 0 ? total / n : 0.0
end

function print_phenotypes(io::IO, g::Generator, input_syms::AbstractVector{Symbol}; pc::PartCatalogue)
    for (i, circuit) in enumerate(iter(g))
        circuit_str = join(circuit, '-')
        output, _ = simulate(Circuit(collect(circuit)), input_syms; pc)
        fun = !any(ismissing, output) ? pheno_label(output) : :oscillator
        println(io, "$i : $circuit_str => $output => $fun")
    end
end

print_phenotypes(g::Generator, input_syms::AbstractVector{Symbol}; pc::PartCatalogue) = print_phenotypes(stdout, g, input_syms; pc)

"""
    characterize(g::Generator, input_syms::AbstractVector{Symbol}; pc::PartCatalogue)

Count the phenotype distribution of all circuits in a generator.
Returns `(pheno_counts, n_circuits)` where `pheno_counts[i]` is the number
of circuits with phenotype `i` (indexed by `pheno_index`).
Any circuits with `missing` outputs are not counted.
"""
function characterize(g::Generator, input_syms::AbstractVector{Symbol}; pc::PartCatalogue)
    k = length(input_syms)
    pheno_counts = zeros(Int, n_phenotypes(k))
    for c in iter(g)
        pheno, _ = simulate(Circuit(collect(c)), input_syms; pc)
        idx = pheno_index(pheno)
        !ismissing(idx) && (pheno_counts[idx] += 1)
    end
    return pheno_counts, number_of_circuits(g)
end

# Get evenly-spaced indices from 1:n
@inline function _get_subsampled_indices(n::Integer, max_sample::Integer)
    n ≤ max_sample && return collect(1:n)
    return round.(Int, range(1, n, length=max_sample))
end

"""
    characterize_with_diversity(g, input_syms; pc, max_sample=50, metric=:edit_distance)

Extends characterize by computing pairwise diversity between circuits per phenotype.

`metric` selects the distance function:
- `:edit_distance` — normalised Levenshtein on functional parts (default)
- `:graph_distance` — active-position normalised Hamming on signed adjacency matrices

Stride-samples up to `max_sample` circuits per phenotype bin.
"""
function characterize_with_diversity(g::Generator, input_syms::AbstractVector{Symbol};
                                     pc::PartCatalogue, max_sample::Int=50,
                                     metric::Symbol=:edit_distance)
    k = length(input_syms)
    n_ph = n_phenotypes(k)
    pheno_counts = zeros(Int, n_ph)

    # Per-phenotype storage: part vectors for edit distance, adjacency matrices for graph distance
    pheno_circuits = [Vector{Symbol}[] for _ in 1:n_ph]
    pheno_matrices = (metric == :graph_distance) ? [Matrix{Int8}[] for _ in 1:n_ph] : nothing

    # Pre-allocate adjacency buffer for graph distance
    adj_buf = (metric == :graph_distance) ? zeros(Int8, length(pc.node_order), length(pc.node_order)) : nothing

    for c in iter(g)
        collected = collect(c)
        parts = collected isa Vector{Symbol} ? collected : collect(Iterators.flatten(collected))
        pheno, _ = simulate(Circuit(parts), input_syms; pc)
        idx = pheno_index(pheno)
        if !ismissing(idx)
            pheno_counts[idx] += 1
            if metric == :graph_distance
                extract_adjacency!(adj_buf, parts; pc)
                push!(pheno_matrices[idx], copy(adj_buf))
            else
                push!(pheno_circuits[idx], parts)
            end
        end
    end
    n_total = number_of_circuits(g)

    # Compute mean pairwise distance per phenotype, then average across phenotypes.
    # Thread across phenotypes because they're independent.
    per_pheno_diversity = fill(-1.0, n_ph)
    Threads.@threads for i in 1:n_ph
        n = (metric == :graph_distance) ? length(pheno_matrices[i]) : length(pheno_circuits[i])

        # Subsample to avoid runtime cost on very large per-phenotype arrays
        idxs = _get_subsampled_indices(n, max_sample)

        # skip unimplemented phenotypes
        length(idxs) < 1 && continue
        # phenotypes with one implementation have diversity 0
        if length(idxs) == 1
            per_pheno_diversity[i] = 0.0
            continue
        end

        pair_sum = 0.0
        n_pairs = 0
        if metric == :graph_distance
            items = pheno_matrices[i]
            for a in idxs, b in idxs
                b > a || continue
                pair_sum += normalised_graph_distance(items[a], items[b])
                n_pairs += 1
            end
        else
            # Strip non-functional parts (i.e. spacers) before computing diversity
            filtered = [filter(∈(pc.functional_parts), pheno_circuits[i][j]) for j in idxs]
            for a in 1:length(filtered), b in a+1:length(filtered)
                pair_sum += normalised_levenshtein(filtered[a], filtered[b])
                n_pairs += 1
            end
        end
        per_pheno_diversity[i] = pair_sum / n_pairs
    end
    # -1 indicates a phenotype is unimplemented
    n_phenos_included = count(≥(0), per_pheno_diversity)
    arch_div = n_phenos_included > 0 ? sum(x for x in per_pheno_diversity if x ≥ 0) / n_phenos_included : 0.0

    return pheno_counts, n_total, arch_div
end


function Base.show(io::IO, g::Generator)
    compact = get(io, :compact, false)
    n_slots = number_of_slots(g)
    n_parts = number_of_parts(g)
    n_circuits = number_of_circuits(g)
    if compact
        show(io, "Generator with $n_parts parts in $(n_slots) slots")
    else
        println(io, "$(typeof(g)) with $n_parts parts in $(n_slots) slots")
        for (i, slot) in enumerate(g.slots)
            print(io, "Slot $i: ")
            !isempty(slot) && print(io, collect(slot))
            println(io)
        end
        println(io, "\nCapable of generating $(n_circuits) circuits")
    end
end

"""
    phenotype_distribution(g::Generator, input_syms; pc::PartCatalogue, show_invalid=true)

Phenotype distribution of `g` as a labelled pair.
Returns `(labels, counts, n_circuits)` in canonical `pheno_labels` order.

`characterize` drops circuits whose simulated output contains `missing` (they
oscillate rather than settling), so its counts do not sum to `n_circuits`. With
`show_invalid=true` those circuits are recovered into a trailing `:invalid`
bucket, making `sum(counts) == n_circuits`.
"""
function phenotype_distribution(g::Generator, input_syms::AbstractVector{Symbol};
                                pc::PartCatalogue, show_invalid::Bool=true)
    k = length(input_syms)
    pheno_counts, n_circuits = characterize(g, input_syms; pc)
    labels = pheno_labels(k)
    counts = pheno_counts
    if show_invalid
        labels = [labels; :invalid]
        counts = [counts; n_circuits - sum(pheno_counts)]
    end
    return labels, counts, n_circuits
end

"""
    summarize(g::Generator, input_syms, pc::PartCatalogue; show_invalid=true)
    summarize(io::IO, g::Generator, input_syms, pc::PartCatalogue; show_invalid=true)

Full summary: structure, building cost, and phenotype distribution barplot.
"""
function summarize(io::IO, g::Generator, input_syms::AbstractVector{Symbol}, pc::PartCatalogue;
                   show_invalid::Bool=true)
    show(io, g)
    println(io, "Building cost: $(get_price(g, pc)) GBP")
    k = length(input_syms)
    labels, counts, n_possible_circuits = phenotype_distribution(g, input_syms; pc)
    n_valid_circuits = n_possible_circuits - counts[end]
    coverage = n_valid_circuits / n_possible_circuits
    println(io)
    println(io, "$n_valid_circuits valid circuits out of $n_possible_circuits")
    println(io, "$n_valid_circuits / $n_possible_circuits = $coverage")
    println(io)

    if !show_invalid
        labels = labels[1:end-1]
        counts = counts[1:end-1]
    end

    # Omit zero-count circuits to reduce barplot size
    if k >= 3
        nonzero = findall(>(0), counts)
        labels = labels[nonzero]
        counts = counts[nonzero]
    end

    println(io, barplot(labels, counts))
end

summarize(g::Generator, input_syms::AbstractVector{Symbol}, pc::PartCatalogue; kwargs...) =
    summarize(stdout, g, input_syms, pc; kwargs...)

# Muted palette; `:invalid` is greyed because it is not a logic function.
const PHENO_FILL   = colorant"#5E7CA0"
const PHENO_EDGE   = colorant"#3F5A75"
const INVALID_FILL = colorant"#B8B8B8"
const INVALID_EDGE = colorant"#9A9A9A"
const PHENO_TEXT   = colorant"#3A3A3A"

"""
    plot_phenotypes(labels, counts; kwargs...)
    plot_phenotypes(g::Generator, input_syms; pc::PartCatalogue, show_invalid=true, kwargs...)

Horizontal bar chart of a phenotype distribution

`kwargs` are forwarded to `Plots.bar`
"""
function plot_phenotypes(labels::AbstractVector, counts::AbstractVector; kwargs...)
    length(labels) == length(counts) || throw(ArgumentError(
        "labels and counts must have equal length, got $(length(labels)) and $(length(counts))"))
    isempty(labels) && throw(ArgumentError("nothing to plot: labels and counts are empty"))

    # Reverse so the first label lands at the top
    ylabels, yvals = string.(reverse(labels)), reverse(counts)
    n = length(ylabels)
    hi = max(maximum(counts), 1)   # an all-zero distribution still needs a finite axis

    defaults = (
        orientation = :h,
        bar_width   = 0.68,
        fillcolor   = [l == "invalid" ? INVALID_FILL : PHENO_FILL for l in ylabels],
        linecolor   = [l == "invalid" ? INVALID_EDGE : PHENO_EDGE for l in ylabels],
        linewidth   = 0.6,
        framestyle  = :axes,
        xshowaxis   = false,
        yticks = (1:n, ylabels),
        ylims  = (0.4, n + 0.6),
        xlims  = (0, hi * 1.16),
        xticks = false,
        xgrid  = false,
        ygrid  = false,
        xlabel = "",
        legend = false,
        tick_direction = :out,
        title = "",
        titlefontsize = 11,
        tickfontsize  = 8,
        size = (700, 500),
        dpi  = 200,
        left_margin   = 4mm,
        right_margin  = 2mm,
        top_margin    = 5mm,
        bottom_margin = 0mm,
        foreground_color_axis   = PHENO_TEXT,
        foreground_color_text   = PHENO_TEXT,
        foreground_color_border = PHENO_TEXT,
    )

    opts = merge(defaults, values(kwargs))
    p = bar(1:n, yvals; opts...)

    # prettier if we redraw the spine over the bars
    xspine = opts.xlims[1]
    plot!(p, [xspine, xspine], collect(opts.ylims);
          color=opts.foreground_color_axis, linewidth=1, label="")

    pad = 0.012 * opts.xlims[2]   # off the drawn axis, so an xlims override still clears the bar
    for (i, v) in enumerate(yvals)
        annotate!(p, v + pad, i, text(string(v), 7, :left, :vcenter, PHENO_TEXT))
    end
    return p
end

function plot_phenotypes(g::Generator, input_syms::AbstractVector{Symbol};
                         pc::PartCatalogue, show_invalid::Bool=true, kwargs...)
    labels, counts, _ = phenotype_distribution(g, input_syms; pc, show_invalid)
    return plot_phenotypes(labels, counts; kwargs...)
end
 