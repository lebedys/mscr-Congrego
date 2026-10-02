#=
Created on 30/01/2023 22:07:41
Last update: -

@author: Michiel Stock
michielfmstock@gmail.com

Make phenotype histograms
=#

using DrWatson
@quickactivate "Congrego"
using Congrego, DataFrames, CSV
using Random, StatsBase, Plots, JSON, UnicodePlots

using Revise
Revise.revise()

# comment this if loaded
#include(scriptsdir("load_precomputed.jl"))

for n_slots in 1:6
    pheno_distr = countmap((logic_functions[v] for (c,v) in gp_map if length(c)==n_slots))
    p = Plots.bar(String.(keys(pheno_distr)), collect(values(pheno_distr)), label="", yscale=:log10)
    Plots.savefig(plotsdir("phenotypes_$(n_slots)_slots.pdf"))
end
