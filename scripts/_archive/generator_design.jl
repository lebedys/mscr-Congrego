#=
Created on 18/01/2023 15:21:03
Last update: 2023-01-19

@author: Michiel Stock
michielfmstock@gmail.com

This script tries to find a good generator construct
that can generate all relevant circuits with a relatively
good distribution
=#

using DrWatson
@quickactivate "Congrego"
using Congrego, Combinatorics, DataFrames, CSV
using Random, StatsBase, Plots, JSON

# load the phenotypes
# -------------------

parts_idict = Dict(p=>i for (i,p) in enumerate(parts))
n_parts = length(parts)

df = CSV.read(datadir("circuits_6_phenotypes.csv"), DataFrame)
# remove cyclic stuff
dropmissing!(df)
# put in a tidy dictionary for querrying
gp_map = Dict(Tuple(Symbol(df[i,j]) for j in 1:6)=>Tuple(df[i,j] for j in 7:10)
                    for i in 1:size(df, 1))

# use matrix to store the phenotypes as binary numbers
# zero means cyclic
const P = zeros(Int, n_parts, n_parts, n_parts, n_parts, n_parts, n_parts);
for (circuit, pheno) in gp_map
    P[(parts_idict[p] for p in circuit)...] = pheno[4] + 2pheno[3] + 4pheno[2] + 8pheno[1] + 1
end

# Generator structure
# -------------------

struct Generator{N,T}
    slots::NTuple{N,Set{T}}
end

Generator(n, parts) = Generator(NTuple{n,Set{eltype(parts)}}([Set([rand(parts)]) for i in 1:n]))
Generator(n=6) = Generator(n, 1:16)


Base.rand(g::Generator) = rand.(g.slots)
iter(g::Generator) = Iterators.product(g.slots...)
Base.collect(g::Generator) = iter(g) |> collect |> vec

number_of_circuits(g::Generator) = prod(length, g.slots)
number_of_slots(g::Generator) = length(g.slots)
number_of_parts(g::Generator) = length(reduce(∪, g.slots))

Base.isless(g1::Generator, g2::Generator) = number_of_circuits(g1) < number_of_circuits(g2)

gen2parts(g::Generator) = [[parts[i] for i in slot] for slot in g.slots]

"""
    characterize(g::Generator)

Performs the phenotypic characterization
of a circuit
"""
function characterize(g::Generator)
    n_circuits = number_of_circuits(g)
    pheno_counts = zeros(Int, 16)
    for c in iter(g)
        ph = P[c...]
        ph > 0 && (pheno_counts[ph] += 1)
    end
    return pheno_counts, n_circuits
end

# GENETIC ALGORITHM
# -----------------

function tournament_selection(objectives, population;
            elitism=true)
    new_pop = empty(population)
    n = length(population)
    while length(new_pop) < n
        i, j = rand(1:n), rand(1:n)
        if objectives[i] > objectives[j]
            push!(new_pop, population[i])
        else
            push!(new_pop, population[j])
        end
    end
    if ellitism
        new_pop[1] = population[argmax(objectives)]
    end
    return new_pop
end

function crossover(g1::Generator, g2::Generator)
    mask = rand(Bool, number_of_slots(g1))
    tup1 = Tuple(mask[i] ? g1.slots[i] : g2.slots[i]
                    for i in 1:number_of_slots(g1))
    tup2 = Tuple(!mask[i] ? g1.slots[i] : g2.slots[i]
                    for i in 1:number_of_slots(g1))
    return Generator(tup1), Generator(tup2)
end

function mutate(g::Generator, parts=1:16; pmut, fdel=sqrt)
    # todo: modify adding/removing prob?
    g = deepcopy(g)
    for slot in g.slots
        rand() ≤ pmut || continue
        # pick a random part to add
        part = rand(parts)
        push!(slot, part)
        l = length(slot)
        # consider deleting a slot
        if l > 1 && rand() ≤ fdel(l / 16)
            delete!(slot, rand(slot))
        else
        end
    end
    return g
end

function recombinate(population; pmut, pcrossover)
    shuffle!(population)
    for i in 1:2:length(population)-1
        if rand() ≤ pcrossover
            population[i:i+1] .= crossover(population[i], population[i+1])
        end
    end
    population .= mutate.(population; pmut)
    return population
end


function genetic_algorithm!(f, population; n_gen,
            pmut=0.05, pcrossover=0.5, elitism=true)
    objectives = zeros(length(population))
    best_obj = zeros(n_gen + 1)
    for t in 1:n_gen
        objectives .= f.(population)
        best_obj[t] = maximum(objectives)
        parents = tournament_selection(objectives, population; elitism)
        population = recombinate(parents; pmut, pcrossover)
    end
    objectives .= f.(population)
    best_obj[end] = maximum(objectives)
    return population, objectives, best_obj
end


"""
Computes entropy, but using a vector
containing only non-oscillating circuits
"""
function funct_entr(g)
    counts, n = characterize(g)
    return StatsBase.entropy(counts ./ n)
end

population = [Generator() for i in 1:100]

population, objectives, best_obj = genetic_algorithm!(funct_entr, population;
                                    n_gen=2000, pcrossover=0.5, pmut=0.025)

plot(best_obj, xlabel="iteration+1", label="funct. entropy",
            title="Generator design GA")
savefig(plotsdir("generator_design_ga.pdf"))

histogram(objectives, xlabel="functi. entropy")

g_best_ga = population[argmax(objectives)]

func_entropy_ga = funct_entr(g_best_ga)
n_circuits_ga = number_of_circuits(g_best_ga)
phenos_g_best_ga, n_circ_ga = characterize(g_best_ga)


open(datadir("generator_ga.json"),"w") do f
    d = Dict(:generator => gen2parts(g_best_ga),
            :phenotypes => phenos_g_best_ga,
            :n_circuits => n_circuits_ga,
            :func_entropy => func_entropy_ga
            )
    JSON.print(f, d)
end

# MAP-elites
# ----------

function descriptors(g::Generator, b=4)
    n_circ = number_of_circuits(g)
    # truncate number any power of b
    n_circ = b^trunc(Int, log(n_circ)/log(b))
    n_parts = number_of_parts(g)
    return n_circ, n_parts
end


function map_elites(f, population; n_gen,
            pmut=0.025, pcrossover=0.5)
    # create archive linking descriptors with solution
    archive = Dict(descriptors(s) => (f(s), s)
                    for s in population)
    best_obj = zeros(n_gen)
    for t in 1:n_gen
        # recombinate the population
        population = recombinate(last.(values(archive)); pmut, pcrossover)
        # store new solutions in archive
        for s in population
            p, b = f(s), descriptors(s)
            if !haskey(archive, b) || p ≥ archive[b][1]
                archive[b] = (p, s)
            end
        end
        best_obj[t] = maximum(first, values(archive))
    end
    return archive, best_obj
end

population = [Generator() for i in 1:100]

map_archive, best_obj = map_elites(funct_entr, population, n_gen=2000)

plot(best_obj, xlabel="iteration", label="funct. entropy",
            title="Generator design MAP-elites")
savefig(plotsdir("generator_design_MAP-elites.pdf"))

# save the stuff
func_entropy_map, g_best_map = maximum(values(map_archive))

n_circuits_map = number_of_circuits(g_best_map)
phenos_g_best_map, n_circ_map = characterize(g_best_map)

map_archive_extended = Dict(d=>(f, gen2parts(g), characterize(g)) for (d, (f, g)) in map_archive)

open(datadir("generator_map-elites.json"), "w") do f
    d = Dict(:generator => gen2parts(g_best_map),
            :phenotypes => phenos_g_best_map,
            :n_circuits => n_circ_map,
            :func_entropy => func_entropy_map,
            :archive => map_archive_extended,
            )
    JSON.print(f, d)
end

trunc_n_circ = map_archive |> keys .|> first |> unique! |> sort!
n_parts_map = map_archive |> keys .|> last |> unique! |> sort!

M = [get(map_archive, (n,p), (missing,missing))[1] for n in trunc_n_circ, p in n_parts_map]
heatmap(M,
    xticks=(1:length(n_parts_map), n_parts_map), xlabel="number of parts used",
    yticks=(1:length(trunc_n_circ), ["[$i, $(4i-1)]" for i in trunc_n_circ]),
    ylabel="number of circuits", title="functional entropy MAP elites")

savefig(plotsdir("generator_design-MAP.png"))
