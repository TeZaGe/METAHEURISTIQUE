using JuMP
using GLPK
using LinearAlgebra
using Random 
using Statistics

include("loadSPP.jl")
include("setSPP.jl")
include("getfname.jl")

# ------------------ Fonction Glouton  Livrable1 ------------------

function utility(ensemble, C, A)
    
    col_cost = sum(A[:, ensemble])
    return col_cost > 0 ? C[ensemble] / col_cost : 0.0
end

function evaluer_solution(solution::Vector{Int}, C)
    if isempty(solution)
        return 0.0
    end
    val = sum(C[i] for i in solution)
    return val
end

function get_lignes_couvertes(solution::Vector{Int,}, A)
    m = size(A, 1) # Nombre de lignes
    lignes_couvertes = zeros(Bool, m)
    for col in solution
        for i in 1:m # On parcourt toutes les lignes
            if A[i, col] == 1
                lignes_couvertes[i] = true
            end
        end
    end
    return lignes_couvertes
end

function check_conflit(candidat::Int, lignes_couvertes::Vector{Bool}, A)
    m = size(A, 1)
    for i in 1:m
        if A[i, candidat] == 1 && lignes_couvertes[i]
            return true # conflit
        end
    end
    return false # Pas de conflit
end

function construction_gloutonne(C,A)
    n = length(C)
    m = size(A, 1)
    solution_final = Int[] 
    solution_utils = Tuple{Int, Float64}[]
    elements = zeros(Bool, m) 

    for i in 1:n
        push!(solution_utils, (i, utility(i, C, A)))
    end
    sort!(solution_utils, by = x -> x[2], rev = true)

    for (index, _) in solution_utils
        # Vérification de conflit 
        a_conflit = false
        for j in 1:m
            if A[j, index] == 1 && elements[j]
                a_conflit = true
                break
            end
        end

        if !a_conflit
            push!(solution_final, index)
            for j in 1:m
                if A[j, index] == 1
                    elements[j] = true
                end
            end
        end
    end
    return solution_final
end

function generer_voisinage_1_1(solution::Vector{Int}, A)
    n = size(A, 2)
    m = size(A, 1)
    
    voisins_set = Set{Vector{Int}}()
    solution_set = Set(solution)
    
    # 1-1 : Échanger un élément avec un candidat
    for s in solution
        
        sol_temp = [x for x in solution if x != s]
        
        #Calculer les lignes couvertes par cette solution de base
        lignes_couvertes_temp = get_lignes_couvertes(sol_temp, A)

        for cand in 1:n
            if cand in solution_set
                continue
            end
            
            if !check_conflit(cand, lignes_couvertes_temp, A)
                # Si pas de conflit, c'est un voisin valide
                new_sol = sort(vcat(sol_temp, cand))
                push!(voisins_set, new_sol)
            end
        end
    end
    
    return collect(voisins_set)
end

function descente_locale(solution_initiale::Vector{Int}, C, A)
    solution_courante = solution_initiale
    amelioration = true
    iteration = 0
    max_iterations = 1000  # Sécurité anti-boucle
    
    while amelioration && iteration < max_iterations
        iteration += 1
        amelioration = false
        
        cout_courant = evaluer_solution(solution_courante, C)
        voisins = generer_voisinage_1_1(solution_courante, A)
        
        meilleur_voisin = nothing
        meilleur_cout = cout_courant
        
        # Évalue TOUS les voisins
        for voisin in voisins
            cout_voisin = evaluer_solution(voisin, C)
            
            if cout_voisin > meilleur_cout
                meilleur_voisin = voisin
                meilleur_cout = cout_voisin
            end
        end
        
        # Si on a trouvé une amélioration
        if meilleur_voisin !== nothing
            solution_courante = meilleur_voisin
            amelioration = true
        end
    end
    
    return solution_courante
end

# ------------------ GRASP Classique et Réactif Livrable2 ------------------

function greedy_randomized_construction(C, 
                                        A, 
                                        all_utilities::Vector{Tuple{Int, Float64}},
                                        α::Float64, 
                                        rcl_size::Int=10)
    
    n = length(C)
    m = size(A, 1)
    elements = zeros(Bool, m)  
    solution_final = Int[]
    available_utilities = copy(all_utilities)
    
    while true
        u_min = Inf
        u_max = -Inf
        valid_candidates = Tuple{Int, Float64}[]

        for (i, util) in available_utilities
            # Vérification de conflit
            a_conflit = false
            for j in 1:m
                if A[j, i] == 1 && elements[j]
                    a_conflit = true
                    break
                end
            end
            
            if !a_conflit
                push!(valid_candidates, (i, util))
                u_min = min(u_min, util)
                u_max = max(u_max, util)
            end
        end
        
        if isempty(valid_candidates); break; end

        uLimit = u_min + α * (u_max - u_min)
        rcl = Int[]
        for (i, util) in valid_candidates
            if util >= uLimit
                push!(rcl, i)
            end
        end

        if isempty(rcl); break; end

        e = rand(rcl) 

        push!(solution_final, e)
        
        for j in 1:m
            if A[j, e] == 1
                elements[j] = true
            end
        end
        
        filter!(x -> x[1] != e, available_utilities)
    end

    return solution_final
end

function select_alpha_index(probabilities::Vector{Float64})
    r = rand()
    cumulative_p = 0.0
    for (i, p) in enumerate(probabilities)
        cumulative_p += p
        if r <= cumulative_p
            return i
        end
    end
    return length(probabilities) 
end

function grasp(A, 
               C, 
               total_iterations::Int,
               alpha::Float64,
               all_utilities::Vector{Tuple{Int, Float64}})
    
    println("--- Démarrage GRASP classique (α = $alpha, total_iterations = $total_iterations) ---")
    
    global_best_sol = Int[]
    global_best_val = -Inf
    
    for iter in 1:total_iterations
        s = greedy_randomized_construction(C, A, all_utilities, alpha, 10)
    
        s_local = descente_locale(s, C, A)
        
        val_local = evaluer_solution(s_local, C)
        
        if val_local > global_best_val
            global_best_val = val_local
            global_best_sol = s_local
            println("It $iter/$total_iterations: Nouv. meilleur global = $global_best_val")
        else
            if iter % 10 == 0
                println("It $iter/$total_iterations: Valeur = $val_local (Meilleur = $global_best_val)")
            end
        end
    end
    
    println("\n--- Fin GRASP classique ---")
    println("Alpha utilisé : $alpha")
    println("Meilleure valeur : $global_best_val")
    
    return global_best_sol, global_best_val
end

function reactive_grasp(A, 
                        C, 
                        total_iterations::Int, 
                        update_block_size::Int,
                        all_utilities::Vector{Tuple{Int, Float64}})
                        
    println("--- Démarrage Reactive GRASP (total_iterations = $total_iterations) ---")

    alphas = [0.1, 0.3, 0.5, 0.7, 0.9]
    num_alphas = length(alphas)
    
    probabilities = ones(num_alphas) ./ num_alphas 
    alpha_scores_sum = zeros(num_alphas) 
    alpha_counts = zeros(Int, num_alphas)     

    global_best_sol = Int[]
    global_best_val = -Inf
    
    for iter in 1:total_iterations
        idx = select_alpha_index(probabilities)
        alpha_selected = alphas[idx]

        s = greedy_randomized_construction(C, A, all_utilities, alpha_selected, 10)
        
        s_local = descente_locale(s, C, A)
        
        val_local = evaluer_solution(s_local, C)

        if val_local > global_best_val
            global_best_val = val_local
            global_best_sol = s_local
            println("It $iter/$total_iterations: Nouv. meilleur global = $global_best_val (avec α = $alpha_selected)")
        else
            if iter % 10 == 0 
                 println("It $iter/$total_iterations: Valeur = $val_local (Meilleur = $global_best_val)")
            end
        end

        alpha_scores_sum[idx] += val_local
        alpha_counts[idx] += 1

        if iter % update_block_size == 0
            avg_scores = zeros(num_alphas)
            for i in 1:num_alphas
                if alpha_counts[i] > 0
                    avg_scores[i] = alpha_scores_sum[i] / alpha_counts[i]
                end
            end
            
            k = 5
            qualities = zeros(num_alphas)
            
            if global_best_val > 0 
                for i in 1:num_alphas
                    if avg_scores[i] > 0
                        qualities[i] = (avg_scores[i] / global_best_val)^k
                    end
                end
            end

            total_quality = sum(qualities)
            if total_quality > 0
                probabilities = qualities ./ total_quality
            else
                probabilities = ones(num_alphas) ./ num_alphas
            end
        end
    end
    
    println("\n--- Probabilités finales des α ---")
    for i in 1:num_alphas
        println("  α = $(alphas[i]): P = $(round(probabilities[i], digits=3))")
    end
    
    # Trouver l'alpha avec la meilleure probabilité finale
    best_alpha_idx = argmax(probabilities)
    best_alpha = alphas[best_alpha_idx]
    
    return global_best_sol, global_best_val, best_alpha
end

# ------------------ Algorithme Génétique  Livrable3 ------------------

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

# ------------------ Appel des fonctions ------------------

function resoudreSPP(fname)
    C, A = loadSPP(fname)
    
    solution_gloutonne = construction_gloutonne(C, A)
    valeur_glout = evaluer_solution(solution_gloutonne, C)
    println("Valeur gloutonne (z) : ", valeur_glout)
    solution_finale = descente_locale(solution_gloutonne, C, A)
    valeur_heur = evaluer_solution(solution_finale, C)

    println("Valeur heuristique (^z) : ", valeur_heur)
    println("Solution finale (heuristique) : ", solution_finale)
end

function resoudreGRASP(fname, alpha::Float64, total_iterations::Int)
    C, A = loadSPP(fname)
    n = length(C)
    all_utilities = sort(
        [(i, utility(i, C, A)) for i in 1:n],
        by = x -> x[2],
        rev = true
    )
    sol_grasp, val_grasp = grasp(A, C, total_iterations, alpha, all_utilities)
    println("Valeur GRASP (z*) : ", val_grasp)
    println("Solution GRASP : ", sol_grasp)
end

function resoudreREACRIVEGRASP(fname, total_iterations::Int, update_block_size::Int)
    C, A = loadSPP(fname)
    n = length(C)
    all_utilities = sort(
        [(i, utility(i, C, A)) for i in 1:n],
        by = x -> x[2],
        rev = true
    )
    sol_reactive, val_reactive, best_alpha = reactive_grasp(A, C, total_iterations, update_block_size, all_utilities)
    println("Valeur Reactive GRASP (z^*) : ", val_reactive)
    println("Solution Reactive GRASP : ", sol_reactive)
    println("Meilleur alpha final : ", best_alpha)
end

function resoudreAG(fname, population_size::Int, generations::Int, mutation_rate::Float64)
    best_individual, best_value = algoGenetique(population_size, generations, mutation_rate, fname)
    println("Valeur AG (z_AG) : ", best_value)
    println("Solution AG : ", best_individual)
end


function experimentationSPP()
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

    for fname in instances
        try
            C, A = loadSPP(fname)
            instance_name = basename(fname)

            # --- Glouton ---
            t_start_glouton = time()
            sol_glouton = construction_gloutonne(C, A)
            t_glouton = time() - t_start_glouton
            z_glouton = evaluer_solution(sol_glouton, C)

            # --- Descente Locale (1-1) ---
            t_start_local = time()
            sol_local = descente_locale(sol_glouton, C, A) 
            t_local = time() - t_start_local
            z_local = evaluer_solution(sol_local, C)

            println(
                "instance = $(instance_name) ",
                "z_glouton = $(round(z_glouton, digits=0)) ",
                "t_glouton = $(round(t_glouton, digits=4)) ",
                "z_local = $(round(z_local, digits=0)) ",
                "t_local = $(round(t_local, digits=4)) ",
            )

        catch e
            println("ERREUR lors du traitement de $fname: $e")
        end
    end
end

function etude_reactivegrasp(; mode="iterations", total_iterations=200, limit_time=60.0)
    """
    - mode="iterations" : chaque run utilise un nombre fixe d'itérations (total_iterations)
    - mode="time" : chaque run s'arrête après limit_time secondes
    """
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
    runs = 3
    update_block = 20
    results = []

    for fname in instances
        println("\nInstance : $fname ")
        C, A = loadSPP(fname)
        n = length(C)
        all_utilities = sort(
            [(i, utility(i, C, A)) for i in 1:n],
            by = x -> x[2],
            rev = true
        )
        
        if mode == "iterations"
            # Mode : nombre fixe d'itérations par run
            println("Mode itérations : $total_iterations itérations par run")
            for run in 1:runs
                t_start = time()
                sol, val, best_alpha = reactive_grasp(A, C, total_iterations, update_block, all_utilities)
                t_end = time()
                cpu_time = t_end - t_start
                
                push!(results, (
                    instance=fname,
                    run=run,
                    best_value=val,
                    cpu_time=cpu_time,
                    alpha_best=best_alpha
                ))
                println("  run = $run | meilleure valeur = $(val) | temps = $(round(cpu_time, digits=2))s | α avec meilleure prob = $best_alpha")
            end
        elseif mode == "time"
            # Mode : limite de temps par run
            println("Mode temps limité : $(limit_time)s maximum par run")
            for run in 1:runs
                t_start = time()
                
                alphas = [0.1, 0.3, 0.5, 0.7, 0.9]
                num_alphas = length(alphas)
                probabilities = ones(num_alphas) ./ num_alphas
                alpha_scores_sum = zeros(num_alphas)
                alpha_counts = zeros(Int, num_alphas)
                global_best_sol = Int[]
                global_best_val = -Inf
                
                iter = 0
                while (time() - t_start) < limit_time
                    iter += 1
                    idx = select_alpha_index(probabilities)
                    alpha_selected = alphas[idx]
                    
                    s = greedy_randomized_construction(C, A, all_utilities, alpha_selected, 10)
                    s_local = descente_locale(s, C, A)
                    val_local = evaluer_solution(s_local, C)
                    
                    if val_local > global_best_val
                        global_best_val = val_local
                        global_best_sol = s_local
                    end
                    
                    alpha_scores_sum[idx] += val_local
                    alpha_counts[idx] += 1
                    
                    if iter % update_block == 0
                        avg_scores = zeros(num_alphas)
                        for i in 1:num_alphas
                            if alpha_counts[i] > 0
                                avg_scores[i] = alpha_scores_sum[i] / alpha_counts[i]
                            end
                        end
                        k = 5
                        qualities = zeros(num_alphas)
                        if global_best_val > 0
                            for i in 1:num_alphas
                                if avg_scores[i] > 0
                                    qualities[i] = (avg_scores[i] / global_best_val)^k
                                end
                            end
                        end
                        total_quality = sum(qualities)
                        if total_quality > 0
                            probabilities = qualities ./ total_quality
                        else
                            probabilities = ones(num_alphas) ./ num_alphas
                        end
                    end
                end
                
                t_end = time()
                cpu_time = t_end - t_start
                sol = global_best_sol
                val = global_best_val
                
                # Trouver l'alpha avec la meilleure probabilité
                best_alpha_idx = argmax(probabilities)
                best_alpha = alphas[best_alpha_idx]
                
                push!(results, (
                    instance=fname,
                    run=run,
                    best_value=val,
                    cpu_time=cpu_time,
                    iterations=iter,
                    alpha_best=best_alpha
                ))
                println("  run = $run | meilleure valeur = $(val) | temps = $(round(cpu_time, digits=2))s | itérations = $iter | α avec meilleure prob = $best_alpha")
            end
        else
            error("Mode inconnu : $mode. Utilisez 'iterations' ou 'time'")
        end
    end

    println("\nRésumé des résultats :")
    if mode == "time"
        println("instance\trun\tbest_value\tcpu_time\titerations\talpha_best")
        for r in results
            println("$(r.instance)\t$(r.run)\t$(r.best_value)\t$(round(r.cpu_time, digits=2))\t$(r.iterations)\t$(r.alpha_best)")
        end
    else
        println("instance\trun\tbest_value\tcpu_time\talpha_best")
        for r in results
            println("$(r.instance)\t$(r.run)\t$(r.best_value)\t$(round(r.cpu_time, digits=2))\t$(r.alpha_best)")
        end
    end
    return results
end

function etude_grasp(; mode="iterations", total_iterations=200, limit_time=60.0)
    """
    Étude de l'influence du paramètre alpha pour GRASP classique.
    - mode="iterations" : chaque run utilise un nombre fixe d'itérations.
    - mode="time" : chaque run s'arrête après un temps limite.
    """
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
    
    alphas_to_test = [0.1, 0.5, 0.9]
    results = []

    for fname in instances
        println("Instance : $fname ")   
        C, A = loadSPP(fname)
        n = length(C)
        all_utilities = sort(
            [(i, utility(i, C, A)) for i in 1:n],
            by = x -> x[2],
            rev = true
        )
        
        for alpha in alphas_to_test
            println("\n--- Test pour alpha = $alpha ---")
            
            if mode == "iterations"
                println("Mode itérations : $total_iterations itérations")
                t_start = time()
                sol, val = grasp(A, C, total_iterations, alpha, all_utilities)
                t_end = time()
                cpu_time = t_end - t_start
                
                push!(results, (
                    instance=fname,
                    alpha=alpha,
                    best_value=val,
                    cpu_time=cpu_time
                ))
                println("  meilleure valeur = $(val) | temps = $(round(cpu_time, digits=2))s")

            elseif mode == "time"
                println("Mode temps limité : $(limit_time)s maximum")
                t_start = time()
                
                global_best_sol = Int[]
                global_best_val = -Inf
                iter = 0
                
                while (time() - t_start) < limit_time
                    iter += 1
                    s = greedy_randomized_construction(C, A, all_utilities, alpha, 10)
                    s_local = descente_locale(s, C, A)
                    val_local = evaluer_solution(s_local, C)
                    
                    if val_local > global_best_val
                        global_best_val = val_local
                        global_best_sol = s_local
                    end
                end
                
                t_end = time()
                cpu_time = t_end - t_start
                
                push!(results, (
                    instance=fname,
                    alpha=alpha,
                    best_value=global_best_val,
                    cpu_time=cpu_time,
                    iterations=iter
                ))
                println("  meilleure valeur = $(global_best_val) | temps = $(round(cpu_time, digits=2))s | itérations = $iter")
            else
                error("Mode inconnu : $mode. Utilisez 'iterations' ou 'time'")
            end
        end
    end

    println("\nRésumé des résultats pour GRASP :")
    if mode == "time"
        println("instance\talpha\tbest_value\tcpu_time\titerations")
        for r in results
            println("$(r.instance)\t$(r.alpha)\t$(r.best_value)\t$(round(r.cpu_time, digits=2))\t$(r.iterations)")
        end
    else
        println("instance\talpha\tbest_value\tcpu_time")
        for r in results
            println("$(r.instance)\t$(r.alpha)\t$(r.best_value)\t$(round(r.cpu_time, digits=2))")
        end
    end
    return results
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

    println("\n### TABLEAU RESULTATS COMPLET (copier-coller) ###")
    println("instance\trun_group\trun_id\trepeat\tpopulation_size\tgenerations\tmutation_rate\tbest_value\tcpu_time")
    for r in all_results
        println("$(r.instance)\t$(r.run_group)\t$(r.run_id)\t$(r.repeat)\t$(r.population_size)\t$(r.generations)\t$(r.mutation_rate)\t$(r.best_value)\t$(round(r.cpu_time, digits=3))")
    end

    return all_results
end

# ----- Instance à résoudre -----
fname = "dat/pb 1000rnd0300.dat"


#------ Appel Livrabel1 ------
resoudreSPP(fname)
# experimentationSPP()

#------ Appel Livrabel2 ------
# resoudreGRASP(fname, 0.4, 200)
# resoudreREACRIVEGRASP(fname, 200, 20)
# etude_grasp(mode="iterations", total_iterations=200)
# etude_reactivegrasp(mode="iterations", total_iterations=200)

#------ Appel Livrabel3 ------
# resoudreAG(fname, 200, 500, 0.4)
# etude_AG(repeats=3)








