#=
Parameter sweep for MAP-Elites generator design.

Usage:
    julia --project=. scripts/sweep_mapelites.jl
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

# ── Fixed constants (not swept) ──────────────────────────────────────────────

const N_SLOTS       = 10
const N_ITER        = 300000
const POP_SIZE      = 200
const COMPOSITE     = true
const PRICE_CAP     = 5000
const MAX_CIRCUITS  = 5000
const MAX_SAMPLE    = 100
const PHENO_WEIGHTS = false
const PCROSSOVER    = 0.25
const PMUT          = 1.0 / N_SLOTS
# distance metric is swept (see sweep_space below)

# Mask of "interesting" 2-input logic functions (logic_functions order), applied
# to the phenotype counts before entropy when PHENO_WEIGHTS is on.
const PHENO_MASK = [
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

# ── Sweep parameter space ────────────────────────────────────────────────────

sweep_space = Dict(
    :w_diversity => [0.0, 0.1, 0.2, 0.3],
    :c_parsimony => [0.0, 0.1, 0.2, 0.3],
    :distance => [:edit_distance, :graph_distance],
)
all_params = dict_list(sweep_space)

# ── Descriptors (shared, no dependency on swept params) ──────────────────────

function descriptors(g::Generator; b=4, price_step=250)
    n_slots = number_of_slots(g)
    n_circ = number_of_circuits(g)
    n_circ = b^trunc(Int, log(n_circ)/log(b))
    n_parts = number_of_parts(g)
    price = get_price(g, pc)
    price = price_step * ceil(Int, price/price_step)
    return n_slots, price
end

# ── Report helpers ────────────────────────────────────────────────────────────

_fmt_dict(d) = join(["$(k) => $(v)" for (k,v) in sort(collect(d), by=x->string(x[1]))], ", ")
_fmt_set(s)  = join(sort(collect(s), by=string), ", ")

# ── Single experiment ─────────────────────────────────────────────────────────

function run_mapelites(swept_params, run_idx, n_total, random_seed)
    # Unpack swept parameters (annotate: dict_list yields Dict{Symbol,Any})
    c_parsimony = swept_params[:c_parsimony]::Float64
    w_diversity = swept_params[:w_diversity]::Float64
    metric      = swept_params[:distance]::Symbol
    w_entropy   = (1 - w_diversity) / 2
    w_coverage  = (1 - w_diversity) / 2
    pmut        = PMUT
    pcrossover  = PCROSSOVER

    Random.seed!(random_seed)
    date = Dates.format(now(), "yyyymmdd_HHMMSS")

    println("  w_diversity=$w_diversity  c_parsimony=$c_parsimony  weights=($w_entropy, $w_coverage, $w_diversity)  pmut=$pmut  seed=$random_seed")

    # ── Filenames ─────────────────────────────────────────────────────────
    dist_tag = metric == :graph_distance ? "graph" : "edit"
    fn_params = (; date, seed = random_seed, dist = dist_tag,
                 wd = w_diversity, wh = w_entropy, wc = w_coverage, p = c_parsimony)
    json_savename = savename("ME", fn_params, "json"; sort=false)
    txt_savename  = savename("ME", fn_params, "txt"; sort=false)
    io = open(datadir("reports", txt_savename), "w")

    # ── Fitness function (captures per-run weights) ───────────────────────
    function fitness(g)
        number_of_circuits(g) > MAX_CIRCUITS && return (0.0, (0.0, 0.0, 0.0, 0.0))
        counts, _, arch_div = characterize_with_diversity(g, external_inputs; pc=pc, max_sample=MAX_SAMPLE, metric=metric)
        n_functions = count(>(0), counts)
        n_functions == 0 && return (0.0, (0.0, 0.0, 0.0, 0.0))

        PHENO_WEIGHTS && (counts[.!PHENO_MASK] .= 0)

        n_inputs = length(external_inputs)
        n_ph = n_phenotypes(n_inputs)
        H = entropy(counts ./ sum(counts))
        raw_e = H / log(n_ph)
        raw_c = log(n_functions) / log(n_ph)
        raw_d = arch_div
        norm_mean_part = clamp(mean_part_bp(g, pc) / MAX_PART_BP, 0.0, 1.0)
        weighted = w_entropy * raw_e + w_coverage * raw_c + w_diversity * raw_d - c_parsimony * norm_mean_part
        return (weighted, (raw_e, raw_c, raw_d, norm_mean_part))
    end

    # ── Population & run ──────────────────────────────────────────────────
    population = COMPOSITE ?
        [CGenerator(N_SLOTS, parts) for _ in 1:POP_SIZE] :
        [Generator(N_SLOTS, parts)  for _ in 1:POP_SIZE]

    timed_result = @timed map_elites(
        fitness,
        population,
        descriptors,
        parts;
        n_iter = N_ITER,
        pmut,
        pcrossover,
        filter_fn = g -> get_price(g, pc) ≤ PRICE_CAP
    )
    map_archive, best_obj, history = timed_result.value

    n_iter_completed = length(best_obj) - 1
    interrupted = n_iter_completed < N_ITER
    interrupted && @warn "Run $run_idx interrupted — captured $n_iter_completed/$N_ITER iterations, $(length(map_archive)) archive cells"

    # ── Progress plot ─────────────────────────────────────────────────────
    plot(
        best_obj,
        xlabel="iteration", label="best fitness",
        title="MAP-Elites sweep run $run_idx/$n_total"
    )
    Plots.savefig(plotsdir(savename("ME_progress", fn_params, "pdf"; sort=false)))

    # ── Metadata ──────────────────────────────────────────────────────────
    run_block = Dict(
        :random_seed           => random_seed,
        :date                  => date,
        :run_index             => run_idx,
        :n_total_runs          => n_total,
        :elapsed_s             => timed_result.time,
        :allocated_bytes       => timed_result.bytes,
        :n_iter_completed      => n_iter_completed,
        :part_catalogue_source => relpath(pc_filepath, projectdir()),
    )
    @tag!(run_block)

    run_metadata = Dict(
        :run => run_block,
        :algorithm => Dict(
            :name                => "map_elites",
            :n_iter              => N_ITER,
            :pop_size            => POP_SIZE,
            :n_slots             => N_SLOTS,
            :pmut                => pmut,
            :pcrossover          => pcrossover,
            :composite           => COMPOSITE,
            :price_cap           => PRICE_CAP,
            :max_circuits        => MAX_CIRCUITS,
            :pheno_weights       => PHENO_WEIGHTS,
            :w_entropy           => w_entropy,
            :w_coverage          => w_coverage,
            :w_diversity         => w_diversity,
            :c_parsimony         => c_parsimony,
            :max_sample          => MAX_SAMPLE,
            :diversity_metric    => metric,
        ),
    )

    # ── Text report ───────────────────────────────────────────────────────
    interrupted && println(io, "\nNOTE: run interrupted after $n_iter_completed/$N_ITER iterations\n")
    println(io, "Sweep run $run_idx/$n_total")
    println(io, "Random seed: $random_seed")
    println(io, "Git: $(run_block[:gitcommit])  script: $(run_block[:script])")
    println(io, "Algorithm: map_elites  n_iter=$N_ITER  pop_size=$POP_SIZE  n_slots=$N_SLOTS  pmut=$pmut  pcrossover=$pcrossover  price_cap=$PRICE_CAP  composite=$COMPOSITE")
    println(io, "Weights: w_entropy=$w_entropy  w_coverage=$w_coverage  w_diversity=$w_diversity")
    println(io, "c_parsimony=$c_parsimony  pheno_weights=$PHENO_WEIGHTS")

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

    _, g_best_map, _ = best_cheapest_generator(map_archive, pc)

    # ── JSON archive ──────────────────────────────────────────────────────
    map_archive_extended = Dict(
        d => let
            (counts, n) = characterize(g, external_inputs; pc=pc)
            Dict(
                :fitness => f,
                :fitness_components => Dict(:entropy => raw_e, :coverage => raw_c, :diversity => raw_d, :norm_mean_part => norm_mean_part),
                :generator => gen2parts(g),
                :phenotypes => counts,
                :n_circuits => n
            )
        end
        for (d, (f, g, (raw_e, raw_c, raw_d, norm_mean_part))) in map_archive if number_of_circuits(g) ≤ MAX_CIRCUITS
    )

    # Serialize history: Dict{descriptor_string => [{iter, fitness, components}, ...]}
    history_json = Dict(
        string(d) => [
            Dict(:iter => t, :fitness => fv, :entropy => raw_e, :coverage => raw_c, :diversity => raw_d, :norm_mean_part => norm_mean_part)
            for (t, fv, (raw_e, raw_c, raw_d, norm_mean_part)) in entries
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

    # ── Text report body ──────────────────────────────────────────────────
    summarize(io, g_best_map, external_inputs, pc)
    print_phenotypes(io, g_best_map, external_inputs; pc=pc)

    println(io, "All other MAP-elite solutions:")
    for k in sort!(collect(keys(map_archive)))
        any(==(0), k) && continue
        f, g, _ = map_archive[k]
        println(io, "\n--- n_slots=$(k[1]), price=$(k[2]) GBP, fitness=$f")
        summarize(io, g, external_inputs, pc)
    end

    # ── Heatmap plots ─────────────────────────────────────────────────────
    n_slots_map = map_archive |> keys .|> first |> unique! |> sort!
    prices = map_archive |> keys .|> last |> unique! |> sort!

    fitness_grid = [get(map_archive, (n, p), (missing, missing, missing))[1] for n in n_slots_map, p in prices]
    Plots.heatmap(fitness_grid,
        xticks=(1:length(prices), prices), xlabel="price of parts (GBP)",
        yticks=(1:length(n_slots_map), n_slots_map),
        ylabel="number of slots used", title="MAP-Elites archive (run $run_idx)")
    Plots.savefig(plotsdir(savename("ME_archive", fn_params, "png"; sort=false)))

    log_circuits_grid = [(n,p) ∈ keys(map_archive) ? log10(number_of_circuits(map_archive[(n,p)][2])) : missing for n in n_slots_map, p in prices]
    Plots.heatmap(log_circuits_grid,
        xticks=(1:length(prices), prices), xlabel="price of parts (GBP)",
        yticks=(1:length(n_slots_map), n_slots_map),
        ylabel="number of slots used", title="log10 circuits (run $run_idx)")
    Plots.savefig(plotsdir(savename("ME_genospace", fn_params, "png"; sort=false)))

    close(io)

    println("  ✓ Completed in $(round(timed_result.time; digits=1))s — archive: $(length(map_archive)) cells, best fitness: $(round(best_obj[end]; digits=4))")
end

# ── Sweep loop ────────────────────────────────────────────────────────────────

Base.exit_on_sigint(false)

for random_seed in [rand(UInt32)]
    println("MAP-Elites parameter sweep")
    println("  $(length(all_params)) runs  (w_diversity ∈ $(sweep_space[:w_diversity]), c_parsimony ∈ $(sweep_space[:c_parsimony]))")
    println("  Random seed: $random_seed")
    println("  Fixed: N_SLOTS=$N_SLOTS  N_ITER=$N_ITER  POP_SIZE=$POP_SIZE  PMUT=$PMUT  PCROSSOVER=$PCROSSOVER  COMPOSITE=$COMPOSITE  PRICE_CAP=$PRICE_CAP  MAX_CIRCUITS=$MAX_CIRCUITS")

    for (i, swept_params) in enumerate(all_params)
        println("\n", "="^60)
        println("Run $i/$(length(all_params)):")
        println("="^60)
        try
            run_mapelites(swept_params, i, length(all_params), random_seed)
        catch e
            if e isa InterruptException
                @warn "Sweep stopped after run $i/$(length(all_params))"
                break
            end
            @error "Run $i failed" exception=(e, catch_backtrace())
        end
    end

    println("\n", "="^60)
    println("Sweep complete.")
    println("="^60)
end
