using JuMP
using GLPK
using LinearAlgebra
using Random 
using Statistics

include("loadSPP.jl")
include("setSPP.jl")
include("getfname.jl")
include("livrableEI1.jl")

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


etude_grasp(mode="time", limit_time=2.0)