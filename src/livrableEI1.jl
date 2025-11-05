using JuMP
using GLPK
using LinearAlgebra
using Random 

include("loadSPP.jl")
include("setSPP.jl")
include("getfname.jl")

function utility(ensemble, C, A)
    col_cost = sum(A[:, ensemble])
    return col_cost > 0 ? C[ensemble] / col_cost : 0.0
end

function evaluer_solution(solution::Vector{Int}, C::Vector{Int})
    if isempty(solution)
        return 0.0
    end
    val = sum(C[i] for i in solution)
    return val
end

function get_lignes_couvertes(solution::Vector{Int}, A::Matrix{Int})
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

function check_conflit(candidat::Int, lignes_couvertes::Vector{Bool}, A::Matrix{Int})
    m = size(A, 1)
    for i in 1:m
        if A[i, candidat] == 1 && lignes_couvertes[i]
            return true # conflit
        end
    end
    return false # Pas de conflit
end


function construction_gloutonne(C::Vector{Int}, A::Matrix{Int})
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

function greedy_randomized_construction(C::Vector{Int}, 
                                        A::Matrix{Int}, 
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
            # Vérification de conflit (simple boucle)
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

function generer_voisinage_1_1(solution::Vector{Int}, A::Matrix{Int})
    n = size(A, 2)
    m = size(A, 1)
    
    voisins_set = Set{Vector{Int}}()
    solution_set = Set(solution)
    
    # 1-1 : Échanger un élément avec un candidat
    for s in solution
        
        # 1. Créer la solution de base (sans 's')
        sol_temp = [x for x in solution if x != s]
        
        # 2. Calculer les lignes couvertes par cette solution de base
        lignes_couvertes_temp = get_lignes_couvertes(sol_temp, A)

        # 3. Essayer d'ajouter chaque candidat (de 1 à n)
        for cand in 1:n
            if cand in solution_set
                continue
            end
            
            # 4. Vérifier le conflit du candidat (simple boucle)
            if !check_conflit(cand, lignes_couvertes_temp, A)
                # Si pas de conflit, c'est un voisin valide
                new_sol = sort(vcat(sol_temp, cand))
                push!(voisins_set, new_sol)
            end
        end
    end
    
    return collect(voisins_set)
end

function descente_locale(solution_initiale::Vector{Int}, C::Vector{Int}, A::Matrix{Int})
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

function reactive_grasp(A::Matrix{Int}, 
                        C::Vector{Int}, 
                        total_iterations::Int, 
                        update_block_size::Int,
                        all_utilities::Vector{Tuple{Int, Float64}})
                        
    println("--- Démarrage Reactive GRASP (total_iterations = $total_iterations) ---")

    alphas = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 0.9]
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
    
    return global_best_sol, global_best_val
end

function resoudreSPP(fname)
    C, A = loadSPP(fname)
    n = length(C)
    
    println("Pré-calcul des Utilités (une seule fois)...")
    all_utilities = sort(
        [(i, utility(i, C, A)) for i in 1:n],
        by = x -> x[2], 
        rev = true
    )

    solution_gloutonne = construction_gloutonne(C, A)
    valeur_glout = evaluer_solution(solution_gloutonne, C)
    println("Valeur gloutonne (z) : ", valeur_glout)

    t_start = time()

    TOTAL_ITERATIONS = 200  
    UPDATE_BLOCK = 20       

    (solution_finale, valeur_heur) = reactive_grasp(
        A, C,
        TOTAL_ITERATIONS,
        UPDATE_BLOCK,
        all_utilities
    )

    t_end = time()
    cpu_time = t_end - t_start

    println("\n--- Fin de resoudreSPP ---")
    println("Valeur heuristique (^z) : ", valeur_heur)
    println("Solution finale (heuristique) : ", solution_finale)
    println("Temps CPU Reactive GRASP (s): ", cpu_time)
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

    # 2. Boucle sur les instances
    for fname in instances
        try
            C, A = loadSPP(fname)
            instance_name = basename(fname)

            # --- 3a. Glouton ---
            t_start_glouton = time()
            sol_glouton = construction_gloutonne(C, A)
            t_glouton = time() - t_start_glouton
            z_glouton = evaluer_solution(sol_glouton, C)

            # --- 3b. Descente Locale (1-1) ---
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
# --- APPEL FINAL ---
# experimentationSPP()


# println("\nLoading...")
fname = "dat/didactic.dat"
solution_heuristique = resoudreSPP(fname)



