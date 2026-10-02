#=
Created on 30/01/2023 17:01:49
Last update: 2023-02-06

@author: Michiel Stock
michielfmstock@gmail.com

Generator design using Genetic Algorithm.
=#

using DrWatson
@quickactivate "Congrego"
using Congrego, DataFrames, CSV
using Random, StatsBase, Plots, JSON, Dates
using UnicodePlots: barplot

# Load part catalogue
include(scriptsdir("michiel_parts.jl"))
# sort for deterministic order
const parts = sort!(collect(pc_m.parts), by=string)
const external_inputs = pc_m.inputs

# see funct_entr below
const pheno_weights = false

n_slots = 6
n_gen = 5000
pop_size = 100
composite = true
pmut = 0.05
target = "entropy_rich"
date = Dates.format(now(), "yyyymmdd_HHMMSS")

seed = rand(UInt32)
Random.seed!(seed)

parameters = (; target, n_gen, pop_size, composite, pheno_weights, date, pmut)

io = open(datadir("reports", savename("design_report_GA", parameters, "txt")), "w")
println(io, "Random seed: $seed")

"""
Computes entropy, but using a vector
containing only non-oscillating circuits
"""
function funct_entr(g; pheno_weights=pheno_weights)
    n_circuits = number_of_circuits(g)
    n_circuits > 2000 && return 0.0
    counts, n = characterize(g, external_inputs; pc=pc_m)
    n_functions = count(>(0), counts)

    # reward only desired functions
    if pheno_weights
        w = [
            false, #:false
            true,  #:and
            false, #:animplyb
            false, #:a
            false, #:bnimplya
            false, #:b
            true,  #:xor
            true,  #:or
            true,  #:nor
            true,  #:xnor
            true,  #:notb
            false, #:bimplya
            true,  #:nota
            false, #:aimplyb
            true,  #:nand
            false, #:true
        ]
        counts[.!w] .= 0
    end 

    obj = exp(0.5entropy(counts ./ n) + 0.5log(n_functions))
    return obj
end

# GENETIC ALGORITHM
# -----------------

println(io,
"""
GENETIC ALGORITHM
-----------------

"""
)

population = composite ? [CGenerator(n_slots, parts) for i in 1:pop_size] : [Generator(n_slots, parts) for i in 1:pop_size];

population, objectives, best_obj = genetic_algorithm!(funct_entr, population, parts;
    n_gen, pcrossover=0.5, pmut)

plot(
    best_obj, 
    xlabel="iteration+1", label="funct. entropy",
    title="Generator design GA"
)
Plots.savefig(plotsdir(savename("GA_progress", parameters, "pdf")))

g_best_ga = population[argmax(objectives)]

func_entropy_ga = funct_entr(g_best_ga)
n_circuits_ga = number_of_circuits(g_best_ga)
phenos_g_best_ga, n_circ_ga = characterize(g_best_ga, external_inputs; pc=pc_m)

open(datadir(savename("generators_GA", parameters, "json")), "w") do f
    d = Dict(:generator => gen2parts(g_best_ga),
        :phenotypes => phenos_g_best_ga,
        :n_circuits => n_circuits_ga,
        :func_entropy => func_entropy_ga,
        :objectives => best_obj,
    )
    JSON.print(f, d)
end

summarize(io, g_best_ga, external_inputs, pc_m)

println(io)
print_phenotypes(io, g_best_ga, external_inputs; pc=pc_m)

close(io)