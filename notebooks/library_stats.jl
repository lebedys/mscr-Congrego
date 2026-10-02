# Phenotype-distribution statistics, shared by ms.ipynb and ms_misassembly.ipynb.

using StatsBase

"""
    library_stats(counts, n_ph)

Size, phenotype coverage and functional entropy of one phenotype distribution.

Coverage and entropy follow `scripts/designer_mapelites.jl`, so they compare with any MAP-Elites
archive; `summarize` prints a different "coverage", valid circuits over possible ones. The
`:invalid` slot is excluded from both and reported on its own.
"""
function library_stats(counts, n_ph)
    named = counts[1:end-1]
    n_functions = count(>(0), named)
    H = entropy(named ./ sum(named))
    return (; size = sum(counts), invalid = counts[end], n_functions,
            coverage = log(n_functions) / log(n_ph),
            entropy = H / log(n_ph), H, perplexity = exp(H))
end

"""
    print_library_stats(name, counts, n_ph; width=20)

One line per library: size, unresolved circuits, coverage and entropy.
"""
function print_library_stats(name, counts, n_ph; width=20)
    s = library_stats(counts, n_ph)
    println(rpad(name, width)
            * lpad(s.size, 6) * " circuits, " * lpad(s.invalid, 5) * " invalid"
            * " ($(rpad(round(100s.invalid / s.size, digits=1), 4))%)"
            * "   coverage $(lpad(s.n_functions, 2))/$n_ph = $(rpad(round(s.coverage, digits=4), 6))"
            * "   entropy $(round(s.entropy, digits=4))"
            * "   ($(round(s.perplexity, digits=2)) effective)")
    return s
end
