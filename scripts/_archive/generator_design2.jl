#=
Created on 25/01/2023 13:58:16
Last update: 2023-01-30

@author: Michiel Stock
michielfmstock@gmail.com

Second version of generator design,
more adaptable so that we can simulate multiple slots.
=#

using DrWatson
@quickactivate "Congrego"
using Congrego, DataFrames, CSV
using Random, StatsBase, Plots, JSON

using Revise
Revise.revise()

include(scriptsdir("load_precomputed.jl"))


# GENETIC ALGORITHM
# -----------------


"""
Computes entropy, but using a vector
containing only non-oscillating circuits
"""
function funct_entr(g)
    n_circuits = number_of_circuits(g)
    n_circuits > 2000 && return 0.0
    counts, n = characterize(g, external_inputs)
    return StatsBase.entropy(counts ./ n)
end

population = [Generator(6) for i in 1:100]

population, objectives, best_obj = genetic_algorithm!(funct_entr, population;
                                    n_gen=2000, pcrossover=0.5, pmut=0.025)

plot(best_obj, xlabel="iteration+1", label="funct. entropy",
            title="Generator design GA")
savefig(plotsdir("generator_design_ga_v2.pdf"))

histogram(objectives, xlabel="funct. entropy")

g_best_ga = population[argmax(objectives)]

func_entropy_ga = funct_entr(g_best_ga)
n_circuits_ga = number_of_circuits(g_best_ga)
phenos_g_best_ga, n_circ_ga = characterize(g_best_ga, external_inputs)


open(datadir("generator_ga_v2.json"),"w") do f
    d = Dict(:generator => gen2parts(g_best_ga),
            :phenotypes => phenos_g_best_ga,
            :n_circuits => n_circuits_ga,
            :func_entropy => func_entropy_ga,
            :objectives => best_obj,
            )
    JSON.print(f, d)
end

# MAP-elites
# ----------

function descriptors(g::Generator; b=4, price_step=250)
    n_slots = number_of_slots(g)
    n_circ = number_of_circuits(g)
    # truncate number any power of b
    n_circ = b^trunc(Int, log(n_circ)/log(b))
    n_parts = number_of_parts(g)
    price = get_price(g)
    price = price_step * ceil(Int, price/price_step)
    return n_slots, price
end

population = [Generator(10) for i in 1:100];

map_archive, best_obj = map_elites(funct_entr, population, descriptors, n_gen=2000);

plot(best_obj, xlabel="iteration", label="funct. entropy",
            title="Generator design MAP-elites")
savefig(plotsdir("generator_design_MAP-elites_v2.pdf"))

# save the stuff
func_entropy_map, g_best_map = maximum(values(map_archive))

n_circuits_map = number_of_circuits(g_best_map)
phenos_g_best_map, n_circ_map = characterize(g_best_map, external_inputs)

map_archive_extended = Dict(d=>(f, gen2parts(g), characterize(g, external_inputs)) for (d, (f, g)) in map_archive
                if number_of_circuits(g) ≤ 5000)

open(datadir("generator_map-elites_v2.json"), "w") do f
    d = Dict(:generator => gen2parts(g_best_map),
            :phenotypes => phenos_g_best_map,
            :n_circuits => n_circ_map,
            :func_entropy => func_entropy_map,
            :archive => map_archive_extended,
            :objectives => best_obj,
            )
    JSON.print(f, d)
end

n_slots_map = map_archive |> keys .|> first |> unique! |> sort!
prices = map_archive |> keys .|> last |> unique! |> sort!

M = [get(map_archive, (n, p), (missing,missing))[1] for n in n_slots_map, p in prices]
heatmap(M,
    xticks=(1:length(prices), prices), xlabel="price of parts (GBP)",
    yticks=(1:length(n_slots_map), n_slots_map),
    ylabel="number of slots used", title="functional entropy MAP elites")

savefig(plotsdir("generator_design-MAP_v2.png"))

Mnc = [haskey(map_archive, (n,p)) ? log10(number_of_circuits(map_archive[(n,p)][2])) : missing for n in n_slots_map, p in prices]
heatmap(Mnc,
    xticks=(1:length(prices), prices), xlabel="price of parts (GBP)",
    yticks=(1:length(n_slots_map), n_slots_map),
    ylabel="number of slots used", title="number of circuits MAP-elites")

savefig(plotsdir("generator_n_circ-MAP.png"))
