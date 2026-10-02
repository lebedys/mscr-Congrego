
# GENETIC ALGORITHM
# -----------------

using Random

export genetic_algorithm!, tournament_selection, crossover,
        mutate, mutate!, recombinate!

"""For a population of n, select winners out of n randomly-chosen pairs of individuals.
If 'elitism' enabled, preserve best individual in first index of return.
"""
function tournament_selection(objectives, population; elitism=true)
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
    if elitism
        new_pop[1] = population[argmax(objectives)]
    end
    return new_pop
end


"Exchange parts in slots between two generators with 50% likelihood."
function crossover(g1::Generator, g2::Generator)
    g1 = deepcopy(g1)
    g2 = deepcopy(g2)
    for (s1, s2) in zip(g1.slots, g2.slots)
        # 50% chance of swapping s1[i] <--> s2[i]
        if !isempty(s1) && !isempty(s2) && rand(Bool)
            p1, p2 = pop!(s1), pop!(s2)
            push!(s1, p2)
            push!(s2, p1)
        end
    end
    return g1, g2
end


mutate(g::Generator, parts; pmut) = mutate(part_type(g), g, parts; pmut)

function mutate(::Type{<:SingletonParts}, g::Generator, parts; pmut)
    # todo: modify adding/removing probabilities?
    g = deepcopy(g)
    for slot in g.slots
        rand() ≤ pmut || continue
        # add / delete / replace with p = 1/3
        if !isempty(slot) && rand() < (1 / 3)
            delete!(slot, rand(slot))
        else
            push!(slot, rand(parts))
            if rand() < (1/2)
                delete!(slot, rand(slot))
            end
        end
    end
    return g
end

function mutate(::Type{<:CompositeParts}, g::Generator, parts; pmut)
    # same as before, but with a twist,
    # a composite part can be made complexer, less complex
    g = deepcopy(g)  # QUESTION: should we just ... use copy?
    for slot in g.slots
        rand() ≤ pmut || continue
        # add / delete / replace with p = 1/3
        if !isempty(slot) && rand() < 1 / 3
            delete!(slot, rand(slot))
        else
            if isempty(slot) || rand(Bool)
                # add a new composite part
                push!(slot, [rand(parts)])
            else
                # modify a part
                comp_part = rand(slot)
                if rand(Bool)
                    # remove a part to the composition
                    if length(comp_part) ≤ 1
                        delete!(slot, comp_part)
                    else
                        pop!(comp_part)
                    end
                else
                    # add a part to the composition
                    push!(comp_part, rand(parts))
                end
            end
        end
    end
    return g
end

mutate!(g::Generator, parts; pmut) = mutate!(part_type(g), g, parts; pmut)

function mutate!(::Type{<:SingletonParts}, g::Generator, parts; pmut)
    for slot in g.slots
        rand() ≤ pmut || continue
        if !isempty(slot) && rand() < (1 / 3)
            delete!(slot, rand(slot))
        else
            push!(slot, rand(parts))
            if rand() < (1/2)
                delete!(slot, rand(slot))
            end
        end
    end
    return g
end

function mutate!(::Type{<:CompositeParts}, g::Generator, parts; pmut)
    for slot in g.slots
        rand() ≤ pmut || continue
        if !isempty(slot) && rand() < 1 / 3
            delete!(slot, rand(slot))
        else
            if isempty(slot) || rand(Bool)
                push!(slot, [rand(parts)])
            else
                comp_part = rand(slot)
                if rand(Bool)
                    if length(comp_part) ≤ 1
                        delete!(slot, comp_part)
                    else
                        pop!(comp_part)
                    end
                else
                    push!(comp_part, rand(parts))
                end
            end
        end
    end
    return g
end

function recombinate!(population::Vector{<:Generator}, parts::AbstractVector{Symbol}; pmut, pcrossover)
    shuffle!(population)
    for i in 1:2:length(population)-1
        if rand() ≤ pcrossover
            population[i:i+1] .= crossover(population[i], population[i+1])
        end
    end
    population .= mutate.(population, Ref(parts); pmut)
    return population
end
    

function genetic_algorithm!(
    f,
    population::Vector{<:Generator},
    parts::AbstractVector{Symbol};
    n_gen, pmut=0.05, pcrossover=0.25, elitism=true, verbose=true)
    t_start = time()
    verbose && println("Starting Genetic Algorithm...")
    objectives = zeros(length(population))
    best_obj = zeros(n_gen + 1)
    for t in 1:n_gen
        Threads.@threads for i in eachindex(population)
            objectives[i] = f(population[i])
        end
        best_obj[t] = maximum(objectives)
        verbose && (t % max(1, n_gen÷100) == 0) &&
            println("gen $t / $n_gen : f_best = $(round(best_obj[t], digits=4)) | elapsed $(round(time()-t_start, digits=1))s")
        parents = tournament_selection(objectives, population; elitism)
        population = recombinate!(parents, parts; pmut, pcrossover)
    end
    Threads.@threads for i in eachindex(population)
        objectives[i] = f(population[i])
    end
    best_obj[end] = maximum(objectives)
    verbose && println("GA complete in $(round(time()-t_start, digits=1))s | f_best = $(round(best_obj[end], digits=4))")
    return population, objectives, best_obj
end

# MAP-elites
# ----------
# TODO adomasl - seems redundant to keep track of fitness, and a vector of its components (e.g. entropy, coverage... etc).
#                Maybe keep only components?

export map_elites

function map_elites(f, population::Vector{<:Generator}, descriptors, parts::AbstractVector{Symbol};
            n_iter, pmut=0.025, pcrossover=0.25, verbose=true, filter_fn=Returns(true))
    t_start = time()
    verbose && println("Starting MAP-Elites...")

    # Infer concrete types for compiler optimisation
    G = eltype(population)
    D = typeof(descriptors(population[1]))

    n_init = length(population)
    # Evaluate first element eagerly to infer concrete types
    init_result = f(population[1])

    # fitness
    init_fits = Vector{Float64}(undef, n_init)
    init_fits[1] = init_result[1]
    
    # individual components of fitness (entropy, coverage, ...etc)
    C = typeof(init_result[2])
    init_fitness_components = Vector{C}(undef, n_init)
    init_fitness_components[1] = init_result[2]
    
    init_descs = Vector{D}(undef, n_init)
    init_descs[1] = descriptors(population[1])

    # Evaluate initial population
    Threads.@threads for i in 2:n_init
        result = f(population[i])
        init_fits[i] = result[1]
        init_fitness_components[i] = result[2]
        init_descs[i] = descriptors(population[i])
    end

    # Seed archive with initial population
    # Note: Dicts are not thread-safe, so update sequentially
    archive = Dict{D, Tuple{Float64, G, C}}()
    for i in 1:n_init
        b = init_descs[i]
        # Note: ties are broken by newer solutions
        if !haskey(archive, b) || init_fits[i] ≥ archive[b][1]
            archive[b] = (init_fits[i], population[i], init_fitness_components[i])
        end
    end

    best_obj = zeros(n_iter + 1)
    best_cell = argmax(v -> v[1], values(archive))
    best_fitness = best_cell[1]
    best_components = best_cell[3]
    best_obj[1] = best_fitness
    last_iter = 0

    # Track fitness history for each archive cell
    history = Dict{D, Vector{Tuple{Int, Float64, C}}}()
    for (b, (f_val, _, c)) in archive
        history[b] = [(0, f_val, c)]
    end

    archive_vals = collect(values(archive))
    # track if archive values change
    archive_updated = false

    try
        for t in 1:n_iter
            verbose && (t % max(1, n_iter÷100) == 0) &&
                println("iter $t / $n_iter : f_best = $(round(best_obj[t], digits=4)) $(round.(best_components, digits=3)) | archive=$(length(archive)) | elapsed $(round(time()-t_start, digits=1))s")

            if archive_updated
                # Rebuild values only when archive changes
                archive_vals = collect(values(archive))
                archive_updated = false
            end

            # Pick random parent from archive
            parent = rand(archive_vals)[2]

            # crossover & mutate
            if rand() ≤ pcrossover
                partner = rand(archive_vals)[2]
                offspring, _ = crossover(parent, partner)
            else
                offspring = deepcopy(parent)
            end
            mutate!(offspring, parts; pmut)

            # Filter out solutions that don't suit our needs, e.g. too costly
            if !filter_fn(offspring)
                best_obj[t + 1] = best_obj[t]
                last_iter = t
                continue
            end

            # Evaluate
            result = f(offspring)
            fitness_val = result[1]
            components_val = result[2]
            descriptor_val = descriptors(offspring)

            # Update archive
            if !haskey(archive, descriptor_val) || fitness_val ≥ archive[descriptor_val][1]
                archive[descriptor_val] = (fitness_val, offspring, components_val)
                archive_updated = true
                if fitness_val ≥ best_fitness
                    best_fitness = fitness_val
                    best_components = components_val
                end
                if !haskey(history, descriptor_val)
                    history[descriptor_val] = Tuple{Int, Float64, C}[]
                end
                push!(history[descriptor_val], (t, fitness_val, components_val))
            end

            best_obj[t + 1] = best_fitness
            last_iter = t
        end
        verbose && println("MAP-Elites complete in $(round(time()-t_start, digits=1))s | archive=$(length(archive)) | f_best = $(round(best_obj[last_iter+1], digits=4))")
    catch e
        # Handle user interrupt
        e isa InterruptException || rethrow()
        @warn "map_elites: interrupted after $last_iter/$n_iter iterations — returning partial archive ($(length(archive)) cells)"
    end
    return archive, best_obj[1:last_iter+1], history
end
