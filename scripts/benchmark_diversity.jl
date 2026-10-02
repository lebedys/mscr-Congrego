#=
Benchmark diversity metrics: :edit_distance vs :graph_distance

Compares the two characterize_with_diversity metrics across generator scales.
Four sections:
  1. extract_adjacency! per circuit (isolated new cost)
  2. Pairwise distance functions head-to-head
  3. characterize_with_diversity end-to-end
  4. max_sample scaling

Usage (run with threads to match production):
    julia --threads=auto --project=. scripts/benchmark_diversity.jl
=#

using Congrego
using BenchmarkTools

include(joinpath(@__DIR__, "parts", "qs2.jl"))
const pc = PC.pc
const input_syms = pc.inputs

n_threads = Threads.nthreads()
println("Julia threads: $n_threads")
n_threads == 1 && @warn "Running single-threaded. characterize_with_diversity uses Threads.@threads — run with `julia --threads=auto` for realistic benchmarks."

# ── Test generators ──────────────────────────────────────────────────────────
# Explicit singleton generators with controlled slot sizes.
# Parts chosen to exercise all _PartClass branches: roadblocker promoters,
# non-RB promoters, TFs, inducible terminators, terminators, spacers.

const generators = [
    ("tiny (2 circuits)", Generator([
        [:pTac, :pLux],
        [:SrpR],
    ])),

    ("small (8 circuits)", Generator([
        [:pTac, :pLux],
        [:SrpR, :BM3R1],
        [:pSrpR, :pBM3R1],
        [:tfSp],
    ])),

    ("medium (72 circuits)", Generator([
        [:pTac, :pLux, :pBetI],
        [:SrpR, :BM3R1],
        [:pSrpR, :pBM3R1, :pLmrA],
        [:BetI, :PhIF],
        [:T10, :TSrpR],
        [:tfSp],
    ])),

    ("large (432 circuits)", Generator([
        [:pTac, :pLux, :pBetI],
        [:SrpR, :BM3R1, :PhIF],
        [:pSrpR, :pBM3R1],
        [:BetI, :LmrA, :AmeR],
        [:T10, :TBM3R1],
        [:pLmrA, :pAmeR],
        [:TSrpR, :tfSp],
        [:tfSp],
    ])),

    ("xlarge (3456 circuits)", Generator([
        [:pTac, :pLux, :pBetI, :pBM3R1],
        [:SrpR, :BM3R1, :PhIF],
        [:pSrpR, :pLmrA, :pAmeR],
        [:BetI, :LmrA, :AmeR],
        [:T10, :TBM3R1],
        [:TSrpR, :tfSp],
        [:pPhIF, :pPsrA],
        [:PhIF, :PsrA],
        [:HlyllR, :AmeR],
        [:tfSp],
    ])),
]

# verify circuit counts
for (label, g) in generators
    n = number_of_circuits(g)
    println("  $label  →  $n circuits ($(number_of_slots(g)) slots, $(number_of_parts(g)) parts)")
end

# ── Section 1: extract_adjacency! per circuit ────────────────────────────────

println("\n", "=" ^ 70)
println("Section 1: extract_adjacency! per circuit")
println("=" ^ 70)

# Representative circuits of different lengths
const test_circuits_adj = [
    ("2-part  [pTac, SrpR]",
     [:pTac, :SrpR]),
    ("4-part  [pTac, SrpR, pSrpR, BM3R1]",
     [:pTac, :SrpR, :pSrpR, :BM3R1]),
    ("6-part  [pTac, SrpR, pSrpR, BM3R1, pBM3R1, tfSp]",
     [:pTac, :SrpR, :pSrpR, :BM3R1, :pBM3R1, :tfSp]),
    ("8-part  [pTac, SrpR, TBM3R1, BetI, pBetI, PhIF, pPhIF, tfSp]",
     [:pTac, :SrpR, :TBM3R1, :BetI, :pBetI, :PhIF, :pPhIF, :tfSp]),
    ("10-part [pTac, SrpR, T10, pSrpR, BM3R1, TBM3R1, BetI, pBetI, PhIF, tfSp]",
     [:pTac, :SrpR, :T10, :pSrpR, :BM3R1, :TBM3R1, :BetI, :pBetI, :PhIF, :tfSp]),
]

N = length(pc.node_order)
adj_buf = zeros(Int8, N, N)

for (label, parts) in test_circuits_adj
    b = @benchmark extract_adjacency!($adj_buf, $parts; pc=$pc)
    m = median(b)
    println("  $(rpad(label, 55))  $(m)")
end

# ── Section 2: pairwise distance functions ───────────────────────────────────

println("\n", "=" ^ 70)
println("Section 2: pairwise distance — graph_distance vs edit_distance")
println("=" ^ 70)

# Build pairs from the medium generator's circuits
println("\nGenerating circuit pairs from medium generator...")
medium_g = generators[3][2]
medium_parts = Vector{Symbol}[]
medium_matrices = Matrix{Int8}[]

for c in iter(medium_g)
    parts = collect(c)
    push!(medium_parts, filter(∈(pc.functional_parts), parts))
    extract_adjacency!(adj_buf, parts; pc=pc)
    push!(medium_matrices, copy(adj_buf))
end

# Pick 10 representative pairs (evenly spaced)
n_circuits = length(medium_parts)
pair_indices = [(i, min(i + div(n_circuits, 10), n_circuits)) for i in 1:div(n_circuits, 10):n_circuits][1:min(10, n_circuits)]

println("  $(length(pair_indices)) pairs from $n_circuits circuits\n")

println("  normalised_levenshtein (per pair):")
for (i, j) in pair_indices
    b = @benchmark normalised_levenshtein($(medium_parts[i]), $(medium_parts[j]))
    println("    pair ($i,$j)  parts=$(length(medium_parts[i]))+$(length(medium_parts[j]))  $(median(b))")
end

println("\n  normalised_graph_distance (per pair):")
for (i, j) in pair_indices
    b = @benchmark normalised_graph_distance($(medium_matrices[i]), $(medium_matrices[j]))
    println("    pair ($i,$j)  $(median(b))")
end

# Aggregate: benchmark a batch of 100 pairwise comparisons
println("\n  Aggregate: 100 pairwise comparisons")
pairs_100 = [(rand(1:n_circuits), rand(1:n_circuits)) for _ in 1:100]

b_lev = @benchmark begin
    s = 0.0
    for (i, j) in $pairs_100
        s += normalised_levenshtein($medium_parts[i], $medium_parts[j])
    end
    s
end
println("    100× normalised_levenshtein:     $(median(b_lev))")

b_gd = @benchmark begin
    s = 0.0
    for (i, j) in $pairs_100
        s += normalised_graph_distance($medium_matrices[i], $medium_matrices[j])
    end
    s
end
println("    100× normalised_graph_distance:  $(median(b_gd))")

# ── Section 3: characterize_with_diversity end-to-end ────────────────────────

println("\n", "=" ^ 70)
println("Section 3: characterize_with_diversity end-to-end")
println("=" ^ 70)
println("  max_sample=50 (default)\n")

println("  $(rpad("Generator", 30))  $(rpad(":edit_distance", 30))  $(rpad(":graph_distance", 30))  ratio")
println("  ", "-"^120)

for (label, g) in generators
    b_ed = @benchmark characterize_with_diversity($g, $input_syms; pc=$pc, metric=:edit_distance) seconds=5 evals=1
    b_gd = @benchmark characterize_with_diversity($g, $input_syms; pc=$pc, metric=:graph_distance) seconds=5 evals=1
    t_ed = median(b_ed).time
    t_gd = median(b_gd).time
    ratio = t_gd / t_ed
    n_samples = min(length(b_ed.times), length(b_gd.times))
    ed_str = "$(round(t_ed / 1e6; digits=2)) ms  ($(median(b_ed).allocs) allocs)"
    gd_str = "$(round(t_gd / 1e6; digits=2)) ms  ($(median(b_gd).allocs) allocs)"
    println("  $(rpad(label, 30))  $(rpad(ed_str, 30))  $(rpad(gd_str, 30))  $(round(ratio; digits=2))×  (n=$n_samples)")
end

# ── Section 4: max_sample scaling ────────────────────────────────────────────

println("\n", "=" ^ 70)
println("Section 4: max_sample scaling (where do the metrics hit a wall?)")
println("=" ^ 70)
println("  Using xlarge generator ($(number_of_circuits(generators[end][2])) circuits)\n")

xlarge_g = generators[end][2]
sample_sizes = [10, 25, 50, 100, 200, 500, 1000, 2000]

println("  $(rpad("max_sample", 14))  $(rpad(":edit_distance", 30))  $(rpad(":graph_distance", 30))  ratio")
println("  ", "-"^110)

for ms in sample_sizes
    b_ed = @benchmark characterize_with_diversity($xlarge_g, $input_syms; pc=$pc, metric=:edit_distance, max_sample=$ms) seconds=10 evals=1
    b_gd = @benchmark characterize_with_diversity($xlarge_g, $input_syms; pc=$pc, metric=:graph_distance, max_sample=$ms) seconds=10 evals=1
    t_ed = median(b_ed).time
    t_gd = median(b_gd).time
    ratio = t_gd / t_ed
    n_samples = min(length(b_ed.times), length(b_gd.times))
    ed_str = "$(round(t_ed / 1e6; digits=1)) ms  ($(median(b_ed).allocs) allocs)"
    gd_str = "$(round(t_gd / 1e6; digits=1)) ms  ($(median(b_gd).allocs) allocs)"
    println("  $(rpad(ms, 14))  $(rpad(ed_str, 30))  $(rpad(gd_str, 30))  $(round(ratio; digits=2))×  (n=$n_samples)")
end

println("\n", "=" ^ 70)
println("Done.")
println("=" ^ 70)
