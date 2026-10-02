#=
Created on 30/01/2023 17:01:49
Last update: 2023-02-06

@author: Michiel Stock
michielfmstock@gmail.com

The master file for designing generators.

update 2023-02-06: add log(richness) to the entropy objective
to maximize the number of circuits
=#

using DrWatson
@quickactivate "Congrego"
using Congrego, DataFrames, CSV
using Random, StatsBase, Plots, JSON, Dates

# Load part catalogue
include(scriptsdir("michiel_parts.jl"))
const external_inputs = pc_m.inputs
const parts = collect(pc_m.parts)

n_gen = 5000
pop_size = 100
composite = true
pheno_weights = false
pmut = 0.05
target = "entropy_rich"
date = today()

pars = (;target, n_gen, pop_size, composite, pheno_weights, date, pmut)

io = open(datadir("reports", savename( "design_report", pars, "txt")), "w")


# GENETIC ALGORITHM
# -----------------

println(io, """
GENETIC ALGORITHM
-----------------

""")

"""
Computes entropy, but using a vector
containing only non-oscillating circuits
"""
function funct_entr(g; pheno_weights=pheno_weights)
    n_circuits = number_of_circuits(g)
    n_circuits > 2000 && return 0.0
    counts, n = characterize(g, external_inputs; pc=pc_m)
    w = [
        false, #Symbol(false), # FFFF
        true, #:and,   # FFFT
        false, #:animplyb,
        false, #:a,
        false, #:bnimplya,
        false, #:b,
        true, #:xor,
        true, #:or,
        true, #:nor,
        true, #:xnor,
        true, #:notb,
        false, #:bimplya,
        true, #nota
        false, #:aimplyb,
        true, #:nand,
        false,#Symbol(true)
    ]
    n_functions = count(>(0), counts)
    pheno_weights && (counts[.!w] .= 0)
    obj = exp(0.5entropy(counts ./ n) + 0.5log(n_functions))
    return obj
end

population = composite ? [CGenerator(6, parts) for i in 1:pop_size] : [Generator(6, parts) for i in 1:pop_size];

population, objectives, best_obj = genetic_algorithm!(funct_entr, population, parts;
                                    n_gen, pcrossover=0.5, pmut)

p=plot(best_obj, xlabel="iteration+1", label="funct. entropy",
            title="Generator design GA")
Plots.savefig(plotsdir(savename("GA_progress", pars, "pdf")))


g_best_ga = population[argmax(objectives)]

func_entropy_ga = funct_entr(g_best_ga)
n_circuits_ga = number_of_circuits(g_best_ga)
phenos_g_best_ga, n_circ_ga = characterize(g_best_ga, external_inputs; pc=pc_m)


open(datadir(savename("generators_GA", pars, "json")),"w") do f
    d = Dict(:generator => gen2parts(g_best_ga),
            :phenotypes => phenos_g_best_ga,
            :n_circuits => n_circuits_ga,
            :func_entropy => func_entropy_ga,
            :objectives => best_obj,
            )
    JSON.print(f, d)
end

summarize(io, g_best_ga, external_inputs, pc_m)

#println(io, barplot(logic_functions, phenos_g_best_ga))

# MAP-elites
# ----------

println(io, """
MAP-ELITES
----------

""")

function descriptors(g::Generator; b=4, price_step=250)
    n_slots = number_of_slots(g)
    n_circ = number_of_circuits(g)
    # truncate number any power of b
    n_circ = b^trunc(Int, log(n_circ)/log(b))
    n_parts = number_of_parts(g)
    price = get_price(g, pc_m)
    price = price_step * ceil(Int, price/price_step)
    return n_slots, price
end

population = composite ?  [CGenerator(10, parts) for i in 1:pop_size] : [Generator(10, parts) for i in 1:pop_size];

map_archive, best_obj = map_elites(funct_entr, 
    population, 
    descriptors, 
    parts;
    n_gen, 
    pmut, 
    filter_fn=g -> get_price(g, pc_m) ≤ 5_000  # avoid very pricy parts
    );

plot(best_obj, xlabel="iteration", label="funct. entropy",
            title="Generator design MAP-elites")
Plots.savefig(plotsdir(savename("MAP-elites_progress", pars, "pdf")))

# save the stuff
func_entropy_map, g_best_map = maximum(values(map_archive))

n_circuits_map = number_of_circuits(g_best_map)
phenos_g_best_map, n_circ_map = characterize(g_best_map, external_inputs; pc=pc_m)

map_archive_extended = Dict(d=>(f, gen2parts(g), characterize(g, external_inputs; pc=pc_m)) for (d, (f, g)) in map_archive
                if number_of_circuits(g) ≤ 5000)

open(datadir(savename("generators_MAP-elites", pars, "json")), "w") do f
    d = Dict(:generator => gen2parts(g_best_map),
            :phenotypes => phenos_g_best_map,
            :n_circuits => n_circ_map,
            :func_entropy => func_entropy_map,
            :archive => map_archive_extended,
            :objectives => best_obj,
            )
    JSON.print(f, d)
end

summarize(io, g_best_map, external_inputs, pc_m)

println(io, "All other MAP-elite solutions:")
for k in sort!(collect(keys(map_archive)))
    println(io, "keys : $k")
    any(==(0), k) && continue
    println(io, map_archive[k])
end

#println(io, barplot(logic_functions, phenos_g_best_map))

n_slots_map = map_archive |> keys .|> first |> unique! |> sort!
prices = map_archive |> keys .|> last |> unique! |> sort!

M = [get(map_archive, (n, p), (missing,missing))[1] for n in n_slots_map, p in prices]
Plots.heatmap(M,
    xticks=(1:length(prices), prices), xlabel="price of parts (GBP)",
    yticks=(1:length(n_slots_map), n_slots_map),
    ylabel="number of slots used", title="effective number of circuits\n MAP elites")

Plots.savefig(plotsdir(savename("MAP-elites_archive", pars, "png")))

Mnc = [haskey(map_archive, (n,p)) ? log10(number_of_circuits(map_archive[(n,p)][2])) : missing for n in n_slots_map, p in prices]
Plots.heatmap(Mnc,
    xticks=(1:length(prices), prices), xlabel="price of parts (GBP)",
    yticks=(1:length(n_slots_map), n_slots_map),
    ylabel="number of slots used", title="log10 number of circuits MAP-elites")

Plots.savefig(plotsdir(savename("MAP-elites_genospace", pars, "png")))

close(io)
