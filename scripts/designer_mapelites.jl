#=
Created on 30/01/2023 17:01:49
Last update: 2023-02-06

@author: Michiel Stock
michielfmstock@gmail.com

Generator design using MAP-Elites.
=#

using DrWatson
@quickactivate "Congrego"
using Congrego, DataFrames, CSV
using Random, StatsBase, Plots, JSON, Dates
using UnicodePlots: barplot

# Load PartCatalogue (pc)
pc_filepath = scriptsdir("parts", "qs2.jl")
include(pc_filepath)
using .PC: pc


# sort for deterministic order
const parts = sort!(collect(pc.parts), by=string)
const external_inputs = pc.inputs

# see fitness below
# TODO adomasl - enable phenotypic weighting (0-100%)
const pheno_weights = false

n_slots    = 10
n_iter     = 200000
pop_size   = 200
composite  = true
pmut       = 1.0/n_slots # 
pcrossover = 0.25
price_cap     = 5000 # skip generators more costly than this (in GBP)
max_circuits  = 5000 # skip generators with more circuits than this

# Fitness function weights
w_diversity = 0.10  #
w_entropy   = (1-w_diversity) / 2  #
w_coverage  = (1-w_diversity) / 2  #
c_parsimony = 0.00  # parsimony penalty coefficient (0 = no penalty)
max_sample  = 100   # per phenotype
metric      = :graph_distance  # :edit_distance or :graph_distance

date = Dates.format(now(), "yyyymmdd_HHMMSS")
seed = rand(UInt32)
Random.seed!(seed)

parameters = (; date, n_iter, w_entropy, w_coverage, w_diversity)
json_savename = savename("mapelites", parameters, "json")
txt_savename = savename("mapelites", parameters, "txt")
io = open(datadir("reports", txt_savename), "w")

"""
Weighted fitness: entropy + coverage + architectural diversity
    
f = w_entropy   * (H / log(n_phenotypes))
  + w_coverage  * (log(n_functions) / log(n_phenotypes))
  + w_diversity * arch_div
"""
function fitness(g; pheno_weights=pheno_weights)
    n_circuits = number_of_circuits(g)
    n_circuits > max_circuits && return (0.0, (0.0, 0.0, 0.0, 0.0))
    counts, _, arch_div = characterize_with_diversity(g, external_inputs; pc=pc, max_sample=max_sample, metric=metric)
    n_functions = count(>(0), counts)
    n_functions == 0 && return (0.0, (0.0, 0.0, 0.0, 0.0))

    k = length(external_inputs)
    n_ph = n_phenotypes(k)

    H = entropy(counts ./ sum(counts))
    raw_e = H / log(n_ph)
    raw_c = log(n_functions) / log(n_ph)
    raw_d = arch_div
    norm_mean_part = clamp(mean_part_bp(g, pc) / MAX_PART_BP, 0.0, 1.0)
    weighted = w_entropy * raw_e + w_coverage * raw_c + w_diversity * raw_d - c_parsimony * norm_mean_part
    return (weighted, (raw_e, raw_c, raw_d, norm_mean_part))
end

function descriptors(g::Generator; b=4, price_step=250)
    n_slots = number_of_slots(g)
    n_circ = number_of_circuits(g)
    # truncate number any power of b
    n_circ = b^trunc(Int, log(n_circ)/log(b))
    n_parts = number_of_parts(g)
    price = get_price(g, pc)
    price = price_step * ceil(Int, price/price_step)
    return n_slots, price
end

population = composite ? [CGenerator(n_slots, parts) for i in 1:pop_size] : [Generator(n_slots, parts) for i in 1:pop_size];

# enable SIGINT to handle early exits with Ctrl-c
Base.exit_on_sigint(false)

# @timed records time & memory allocations
timed_result = @timed map_elites(
    fitness,
    population,
    descriptors,
    parts;
    n_iter,
    pmut,
    pcrossover,
    filter_fn = g -> get_price(g, pc) ≤ price_cap
);
map_archive, best_obj, history = timed_result.value

# Handle user interrupt (Ctrl-c)
n_iter_completed = length(best_obj) - 1
interrupted = n_iter_completed < n_iter
interrupted && @warn "Run interrupted — captured $n_iter_completed/$n_iter iterations, $(length(map_archive)) archive cells"

plot(
    best_obj,
    xlabel="iteration", label="funct. entropy",
    title="Generator design MAP-Elites"
)
Plots.savefig(plotsdir(savename("MAP-Elites_progress", parameters, "pdf")))

# Build run metadata (post-run: requires timed_result)
run_block = Dict(
    :random_seed           => seed,
    :date                  => date,
    :elapsed_s             => timed_result.time,
    :allocated_bytes       => timed_result.bytes,
    :n_iter_completed      => n_iter_completed,
    :part_catalogue_source => relpath(pc_filepath, projectdir()),
)
@tag!(run_block)   # adds :gitcommit, :gitpatch, :script

run_metadata = Dict(
    :run => run_block,
    :algorithm => Dict(
        :name                => "map_elites",
        :n_iter              => n_iter,
        :pop_size            => pop_size,
        :n_slots             => n_slots,
        :pmut                => pmut,
        :pcrossover          => pcrossover,
        :composite           => composite,
        :price_cap           => price_cap,
        :archive_circuit_cap => 5000,  # archive entries with n_circuits > this are excluded
        :pheno_weights       => pheno_weights,
        :w_entropy           => w_entropy,
        :w_coverage          => w_coverage,
        :w_diversity         => w_diversity,
        :c_parsimony         => c_parsimony,
        :max_sample          => max_sample,
        :diversity_metric    => metric,
    ),
)
# input_symbols intentionally NOT stored — use part_catalogue.inputs as single source of truth

# Print metadata header to txt report
interrupted && println(io, "\nNOTE: run interrupted after $n_iter_completed/$n_iter iterations\n")
println(io, "Random seed: $seed")
println(io, "Git: $(run_block[:gitcommit])  script: $(run_block[:script])")
println(io, "Algorithm: map_elites  n_iter=$n_iter  pop_size=$pop_size  n_slots=$n_slots  pmut=$pmut  pcrossover=$pcrossover  price_cap=$price_cap  composite=$composite")
println(io, "pheno_weights=$pheno_weights")

# Full part catalogue — one line per field
_fmt_dict(d) = join(["$(k) => $(v)" for (k,v) in sort(collect(d), by=x->string(x[1]))], ", ")
_fmt_set(s)  = join(sort(collect(s), by=string), ", ")
println(io, "\nPart catalogue: $(run_block[:part_catalogue_source])  ($(length(pc.parts)) parts)")
println(io, "  parts                 = { $(_fmt_set(pc.parts)) }")
println(io, "  inputs                = $(pc.inputs)")
println(io, "  activatable_promoters = { $(_fmt_dict(pc.activatable_promoters)) }")
println(io, "  repressible_promoters = { $(_fmt_dict(pc.repressible_promoters)) }")
println(io, "  terminators           = { $(_fmt_set(pc.terminators)) }")
println(io, "  roadblockers          = { $(_fmt_set(pc.roadblockers)) }")
println(io, "  inducible_terminators = { $(_fmt_dict(pc.inducible_terminators)) }")
println(io, "  transcription_factors = { $(_fmt_set(pc.transcription_factors)) }")

println(io,
"""
MAP-ELITES
----------

"""
)

func_entropy_map, g_best_map, _ = best_cheapest_generator(map_archive, pc)

n_circuits_map = number_of_circuits(g_best_map)
phenos_g_best_map, n_circ_map = characterize(g_best_map, external_inputs; pc=pc)

# Archive entries with named fields
map_archive_extended = Dict(
    d => let
        (counts, n) = characterize(g, external_inputs; pc=pc)

        Dict(
            :fitness => f,
            :fitness_components => Dict(:entropy => c[1], :coverage => c[2], :diversity => c[3], :norm_mean_part => c[4]),
            :generator => gen2parts(g),
            :phenotypes => counts,
            :n_circuits => n
        )
        end
    for (d, (f, g, c)) in map_archive if number_of_circuits(g) ≤ max_circuits
)

# Serialize history: Dict{descriptor_string => [{gen, fitness, components}, ...]}
history_json = Dict(
    string(d) => [
        Dict(:iter => t, :fitness => fv, :entropy => c[1], :coverage => c[2], :diversity => c[3], :norm_mean_part => c[4])
        for (t, fv, c) in entries
    ]
    for (d, entries) in history
)

open(datadir(json_savename), "w") do f
    JSON.print(f, Dict(
        :archive                    => map_archive_extended,
        :part_catalogue             => pc,
        :max_fitness_per_iteration  => best_obj,
        :archive_history            => history_json,
        :metadata                   => run_metadata,
    ))
end

summarize(io, g_best_map, external_inputs, pc)
print_phenotypes(io, g_best_map, external_inputs; pc=pc)

println(io, "All other MAP-elite solutions:")
for k in sort!(collect(keys(map_archive)))
    any(==(0), k) && continue
    f, g, _ = map_archive[k]
    println(io, "\n--- n_slots=$(k[1]), price=$(k[2]) GBP, entropy=$f")
    summarize(io, g, external_inputs, pc)
end

#println(io, barplot(logic_functions, phenos_g_best_map))

n_slots_map = map_archive |> keys .|> first |> unique! |> sort!
prices = map_archive |> keys .|> last |> unique! |> sort!

M = [get(map_archive, (n, p), (missing, missing, missing))[1] for n in n_slots_map, p in prices]
Plots.heatmap(M,
    xticks=(1:length(prices), prices), xlabel="price of parts (GBP)",
    yticks=(1:length(n_slots_map), n_slots_map),
    ylabel="number of slots used", title="effective number of circuits\n MAP elites")

Plots.savefig(plotsdir(savename("MAP-elites_archive", parameters, "png")))

Mnc = [haskey(map_archive, (n,p)) ? log10(number_of_circuits(map_archive[(n,p)][2])) : missing for n in n_slots_map, p in prices]
Plots.heatmap(Mnc,
    xticks=(1:length(prices), prices), xlabel="price of parts (GBP)",
    yticks=(1:length(n_slots_map), n_slots_map),
    ylabel="number of slots used", title="log10 number of circuits MAP-elites")

Plots.savefig(plotsdir(savename("MAP-elites_genospace", parameters, "png")))

close(io)
