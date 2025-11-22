using JuMP
using GLPK
using LinearAlgebra
using Random 
using Statistics

include("loadSPP.jl")
include("setSPP.jl")
include("getfname.jl")
include("livrableEI1.jl")

function selectionParent(population::Vector{Vector{Int}}, C)
    candidates = rand(1:length(population), 2)
    best_idx = candidates[1]
    best_val = evaluer_solution(population[best_idx], C)
    for idx in candidates[2:end]
        v = evaluer_solution(population[idx], C)
        if v > best_val
            best_idx = idx
            best_val = v
        end
    end
    return population[best_idx]
end

function crossover_set(parent1::Vector{Int}, parent2::Vector{Int}, A)
    s1 = Set(parent1); s2 = Set(parent2) 
    child = collect(intersect(s1,s2)) 
    ligne_couverte = get_lignes_couvertes(child, A) 

    # Ajout aléatoire d'éléments des deux parents
    candidates = shuffle(collect(setdiff(union(s1,s2), Set(child))))

    for c in candidates
        if !check_conflit(c, ligne_couverte, A) 
            push!(child, c)

            for i in 1:size(A,1)
                if A[i,c] == 1
                    ligne_couverte[i] = true 
                end
            end
        end
    end
    sort!(child)
    return child
end

function mutation_set(individual::Vector{Int}, n::Int, A, mutation_rate::Float64)
    if rand() < mutation_rate
        ligne_couverte_lines = get_lignes_couvertes(individual, A)
        if rand() < 0.5 && !isempty(individual)
            popat!(individual, rand(1:length(individual)))
        else
            cand = rand(1:n)
            
            if !(cand in individual) 
                if !check_conflit(cand, ligne_couverte_lines, A)
                    push!(individual, cand)
                end
            end
        end
        sort!(individual)
    end
    return individual
end

function survivantEnfant(enfant1, enfant2, C)
    val1 = evaluer_solution(enfant1, C)
    val2 = evaluer_solution(enfant2, C)

    if val1 > val2
        return enfant1
    else
        return enfant2
    end

end

function algoGenetique(population_size, generations, mutation_rate, fname)
    C, A = loadSPP(fname)
    n = length(C)

    all_utilities = sort(
        [(i, utility(i, C, A)) for i in 1:n],
        by = x -> x[2], 
        rev = true
    )

    println("Init population")
    population = Vector{Vector{Int}}(undef, population_size)
    for i in 1:population_size
        alpha_aleatoire = 0.3 + rand() * 0.5 
        population[i] = greedy_randomized_construction(C, A, all_utilities, alpha_aleatoire, 10)
    end

    best_individual = population[1]
    best_value = evaluer_solution(best_individual, C)
    
    for i in 2:population_size 
        individu_actuel = population[i] 
        val = evaluer_solution(individu_actuel, C)
        
        if val > best_value
            best_value = val
            best_individual = individu_actuel
        end
    end
    println("Meilleur score initial : $best_value")

    # --- BOUCLE D'ÉVOLUTION ---
    for gen in 1:generations
        new_population = Vector{Vector{Int}}()

        #ÉLITISME
        push!(new_population, best_individual)

        while length(new_population) < population_size
            parent1 = selectionParent(population, C; )
            parent2 = selectionParent(population, C;)

            enfant1 = crossover_set(parent1, parent2, A)
            enfant2 = crossover_set(parent2, parent1, A)

            enfant1 = mutation_set(enfant1, n, A, mutation_rate)
            enfant2 = mutation_set(enfant2, n, A, mutation_rate)

            survivant = survivantEnfant(enfant1, enfant2, C)
            push!(new_population, survivant)

            val_survivant = evaluer_solution(survivant, C)
            if val_survivant > best_value
                best_value = val_survivant
                best_individual = survivant
                if gen % 10 == 0 # Affiche tous les 10 générations
                    println("Gen $gen: Nouveau meilleur score global = $best_value")
                end
            end
        end
        population = new_population
    end

    println("Fin de l'AG. Meilleur score final : $best_value")
    return best_individual, best_value
end


function etude_AG(; repeats::Int)
    instances = [
        "dat/pb_100rnd0100.dat",
        "dat/pb_200rnd0100.dat",
        "dat/pb_500rnd0300.dat",
        "dat/pb_500rnd1500.dat",
        "dat/pb_500rnd1700.dat",
        "dat/pb_200rnd0400.dat",
        "dat/pb_200rnd0700.dat",
        "dat/pb_1000rnd0100.dat",
        "dat/pb_1000rnd0300.dat",
        "dat/didactic.dat"
    ]

    base_pop = 200
    base_gen = 500
    base_mut = 0.5

    pop_sizes = [100, 150, 200]     
    gens      = [100, 250, 500]         
    muts      = [0.1, 0.5, 0.7]         

    all_results = []

    println("Starting etude_AG with repeats=$repeats")
    for fname in instances
        println("\n=== Instance : $fname ===")
        results = []
        run_idx = 1

        # On change la population_size
        for ps in pop_sizes
            for rep in 1:repeats
                t_start = time()
                _, val = algoGenetique(ps, base_gen, base_mut, fname)
                t_end = time()
                push!(results, (
                    instance = fname,
                    run_group = "pop_size",
                    run_id = run_idx,
                    repeat = rep,
                    population_size = ps,
                    generations = base_gen,
                    mutation_rate = base_mut,
                    best_value = val,
                    cpu_time = t_end - t_start
                ))
            end
            run_idx += 1
        end

        # on change les generations size
        for g in gens
            for rep in 1:repeats
                t_start = time()
                _, val = algoGenetique(base_pop, g, base_mut, fname)
                t_end = time()
                push!(results, (
                    instance = fname,
                    run_group = "generations",
                    run_id = run_idx,
                    repeat = rep,
                    population_size = base_pop,
                    generations = g,
                    mutation_rate = base_mut,
                    best_value = val,
                    cpu_time = t_end - t_start
                ))
            end
            run_idx += 1
        end

        # on change les mutation rates
        for mu in muts
            for rep in 1:repeats
                t_start = time()
                _, val = algoGenetique(base_pop, base_gen, mu, fname)
                t_end = time()
                push!(results, (
                    instance = fname,
                    run_group = "mutation",
                    run_id = run_idx,
                    repeat = rep,
                    population_size = base_pop,
                    generations = base_gen,
                    mutation_rate = mu,
                    best_value = val,
                    cpu_time = t_end - t_start
                ))
            end
            run_idx += 1
        end
        append!(all_results, results)
    end

    println("instance\trun_group\trun_id\trepeat\tpopulation_size\tgenerations\tmutation_rate\tbest_value\tcpu_time")
    for r in all_results
        println("$(r.instance)\t$(r.run_group)\t$(r.run_id)\t$(r.repeat)\t$(r.population_size)\t$(r.generations)\t$(r.mutation_rate)\t$(r.best_value)\t$(round(r.cpu_time, digits=3))")
    end

    return all_results
end


algoGenetique(200, 500, 0.5, "dat/pb_500rnd1500.dat")