#=
Compare Congrego's simulate() against Sim's precomputed phenotypes.
=#

using DrWatson
@quickactivate "Congrego"
using Congrego, DataFrames, CSV, JSON

include(scriptsdir("michiel_parts.jl"))


#= This CSV contains all of the phenotypes in the library from Sim's thesis, modelled by his original simulator.
All circuits have 6 parts.

-- INPUT ORDERING --
Input columns are ordered as follows:
  |            -/-            |            -/+           |           +/-            |           +/+           |
  |   IPTG=false, aTC=false   |   IPTG=false, aTC=true   |   IPTG=true, aTC=false   |   IPTG=true, aTC=true   |


This ordering can be deduced from the circuit on line 17, where we see the states for the two promoters induced accordingly:
    17,pTet,pBAD,AraC,T10,pTac,tfSp,"[[0], [0], [0], [0], [0], [0]]","[[1], [0], [0], [0], [0], [0]]","[[0], [0], [0], [0], [1], [1]]","[[1], [0], [0], [0], [1], [1]]","[[0], [0], [1], [1]]"

    Circuit: pTet-pBAD-AraC-T10-pTac-tfSp
        |            -/-            |            -/+           |           +/-            |           +/+           |
        |   IPTG=false, aTC=false   |   IPTG=false, aTC=true   |   IPTG=true, aTC=false   |   IPTG=true, aTC=true   |
  pTet  |           false           |            true          |          false           |           true          |
  pTac  |           false           |           false          |           true           |           true          |


-- EXAMPLES OF OTHER CIRCUITS --
```
,part0,part1,part2,part2,part4,part5,-/-,-/+,+/-,+/+,output_phenotype
0,pTet,pBAD,AraC,pTet,pBAD,AraC,"[[0], [0], [0], [0], [0], [0]]","[[1], [0], [0], [1], [0], [0]]","[[0], [0], [0], [0], [0], [0]]","[[1], [0], [0], [1], [0], [0]]","[[0], [0], [0], [0]]"
...
30,pTet,pBAD,PsrA,T10,pPsrA,AraC,"[[0], [0], [0], [0], [1], [1]]","[[1], [0], [0], [0], [1], [1]]","[[0], [0], [0], [0], [1], [1]]","[[1], [0], [0], [0], [1], [1]]","[[1], [1], [1], [1]]"
...
319,T10,pTac,tfSp,T10,pPsrA,PsrA,"[[0, 0], [0, 0], [0, 0], [0, 0], [1, 0], [1, 0]]","[[0, 0], [0, 0], [0, 0], [0, 0], [0, 1], [0, 1]]","[[0, 0], [1, 1], [1, 1], [0, 0], [1, 0], [1, 0]]","[[0, 0], [1, 1], [1, 1], [0, 0], [0, 1], [0, 1]]","[[1, 0], [0, 1], [1, 0], [0, 1]]"
...
```
=#
const reference_csv = joinpath("..", "generated_data_ms", "data", "library_full_df-sim.csv")

function parse_phenotype(s)
    parsed = JSON.parse(s)
    ntuple(4) do i
        vals = parsed[i]
        val  = vals[1]
        if length(vals) > 1 && !all(==(val), vals)
            missing   # cycling gfp values → oscillation
        else
            isnothing(val) ? missing : Bool(val)
        end
    end
end

references = CSV.read(reference_csv, DataFrame)
part_cols  = names(references)[2:7]   # part0–part5; positional to avoid duplicate-name issue
diff_dict  = Dict{String, Vector{String}}()

#=
Formatted example:
    "pTet-pPsrA-AraC-T10-pBAD-PsrA": [
        "0 0 0 0",   <--- "expected" (Sim's)
        "- 1 - 1"    <--- "got"      (Congrego)
    ],

Remember that the input order is:
|   IPTG=false, aTC=false   |   IPTG=false, aTC=true   |   IPTG=true, aTC=false   |   IPTG=true, aTC=true   |
=#
fmt(tup) = join([ismissing(v) ? "-" : (v ? "1" : "0") for v in tup], " ")
for row in eachrow(references)
    parts    = Symbol.([row[col] for col in part_cols])
    expected = parse_phenotype(row[:output_phenotype])
    got, _   = simulate(Circuit(parts), [:IPTG, :aTC]; pc=pc_m)
    if !isequal(collect(got), collect(expected))
        diff_dict[join(parts, "-")] = [fmt(expected), fmt(got)]
    end
end

println("$(nrow(references)) circuits checked, $(length(diff_dict)) differences")

open(datadir("Congrego-Sim-diff.json"), "w") do f
    JSON.print(f, diff_dict)
end
