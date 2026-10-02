#=
Benchmark simulate()
=#

using Congrego
using BenchmarkTools

include(joinpath(@__DIR__, "parts", "ms_original.jl"))
const pc = PC.pc

# test circuits: short to long
const circuits = [
    Circuit(Symbol.(split("pTac-tfSp", "-"))),                                          # 2 parts
    Circuit(Symbol.(split("pTac-SrpR-pSrpR-BM3R1", "-"))),                              # 4 parts
    Circuit(Symbol.(split("pTet-PhIF-pPhIF-SrpR-pSrpR-BM3R1", "-"))),                   # 6 parts
    Circuit(Symbol.(split("pTac-SrpR-pSrpR-BM3R1-pBM3R1-PhIF-pPhIF-tfSp", "-"))),       # 8 parts
    Circuit(Symbol.(split("pTac-SrpR-T10-pTet-PhIF-pPhIF-BM3R1-pBM3R1-LmrA-pLmrA-tfSp", "-"))), # 11 parts
    Circuit(Symbol.(split("pTet-pTac-TSrpR-BM3R1-LmrA-AraC-PhIF-pSrpR-pTra-BetI-pTet-BM3R1-IcaRA-pBM3R1-pBAD", "-"))), # 15 parts
    Circuit(Symbol.(split("pTac-SrpR-pSrpR-BM3R1-pBM3R1-PhIF-pPhIF-LmrA-pLmrA-AmeR-pAmeR-tfSp-T10-pTet-IcaRA-pIcaRA-BetI-pBetI-tfSp", "-"))), # 19 parts
    Circuit(Symbol.(split("pTac-SrpR-pSrpR-BM3R1-pBM3R1-PhIF-pPhIF-LmrA-pLmrA-AmeR-pAmeR-tfSp-T10-pTet-IcaRA-pIcaRA-BetI-pBetI-QacR-pQacR-LitR-pLitR-tfSp", "-"))), # 23 parts
]

const input_syms = [:IPTG, :aTC]

println("=" ^ 60)
println("simulate() benchmark")
println("=" ^ 60)

# single-combo benchmarks
println("\n--- single input combo: simulate(circuit, :IPTG => true, :aTC => false) ---")
for c in circuits
    label = join(string.(c.parts), "-")
    b = @benchmark simulate($c, :IPTG => true, :aTC => false; pc=$pc)
    println("  $(rpad(label, 55))  $(median(b))")
end

# all-combos benchmarks
println("\n--- all combos: simulate(circuit, [:IPTG, :aTC]) ---")
for c in circuits
    label = join(string.(c.parts), "-")
    b = @benchmark simulate($c, $input_syms; pc=$pc)
    println("  $(rpad(label, 55))  $(median(b))")
end

# bulk benchmark: simulate many circuits (typical characterize workload)
println("\n--- bulk: 100 circuits × all combos ---")
const bulk_circuits = [circuits[rand(1:end)] for _ in 1:100]
b = @benchmark begin
    for c in $bulk_circuits
        simulate(c, $input_syms; pc=$pc)
    end
end
println("  100 circuits × 4 combos:  $(median(b))")

println("\n--- allocation check: single call ---")
c = circuits[4]
label = join(string.(c.parts), "-")
b = @benchmark simulate($c, :IPTG => true, :aTC => false; pc=$pc)
println("  $label")
println("  median time:   $(median(b).time / 1e3) μs")
println("  allocations:   $(median(b).allocs)")
println("  memory:        $(median(b).memory) bytes")
