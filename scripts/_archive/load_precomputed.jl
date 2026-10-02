#=
Created on 27/01/2023 11:23:58
Last update: 2023-02-01

@author: Michiel Stock
michielfmstock@gmail.com

Load the precomputed cicuits.

Update: as of 1/2/2023 this is deprecated. The space is too big to simulate them all anyway
=#

using Congrego, DataFrames, CSV

# load the phenotypes
# -------------------

parts_idict = Dict(p=>i for (i,p) in enumerate(parts))
n_parts = length(parts)

# put in a tidy dictionary for querrying
const gp_map = Dict{Vector{Int},Int}(Int[]=>1)
for n_slots in 1:6

    df = CSV.read(datadir("circuits_$(n_slots)_phenotypes.csv"), DataFrame)
    # remove cyclic stuff
    dropmissing!(df)

    for i in 1:nrow(df)
        circuit = [parts_idict[Symbol(df[i,j])] for j in 1:n_slots]
        pheno = [df[i,j] for j in range(n_slots+1, length=4)]
        gp_map[circuit] = pheno[4] + 2pheno[3] + 4pheno[2] + 8pheno[1] + 1
    end
end
