
using JuMP
using GLPK
using LinearAlgebra

# using CSV
# using DataFrames
include("loadSPP.jl")
include("setSPP.jl")

include("getfname.jl")


#pb_200rnd0700
#pb_500rnd1500
#pb_1000rnd0300 z_glouton = 507 ts = 0,91 z_heuristique = 538 ts = 2,59
#pb_2000rnd0700
println("\nLoading...")
fname = "Data/pb_1000rnd0300.dat"
C, A = loadSPP(fname)



# le ! dans les fonction pour indiquer que la fonction modifie ses arguments

function resoudreSPP(fname)
    C, A = loadSPP(fname)
    # @show C
    # @show A

    
    t_start = time()

    # heuristique de construction gloutonne
    solution_initiale = construction_gloutonne(C, A)
    println("Solution initiale (gloutonne) : ", solution_initiale)

    valeur_glout = sum(C[i] for i in solution_initiale)
    println("Valeur gloutonne (z) : ", valeur_glout)

    # amélioration par descente locale
    solution_finale = descente_locale(solution_initiale, C, A)
    valeur_heur = sum(C[i] for i in solution_finale)

    t_end = time()
    cpu_time = t_end - t_start

    println("Valeur heuristique (^z) : ", valeur_heur)
    println("Solution finale (heuristique) : ", solution_finale)
    println("Temps CPU total (s): ", cpu_time)


    # time_start = time()

    # grasp_solution = grasp(A, C, 5)
    # grasp_value = sum(C[i] for i in grasp_solution)
    # time_end = time()
    # grasp_cpu_time = time_end - time_start
    # println("Valeur GRASP (z) : ", grasp_value)
    # println("Solution GRASP : ", grasp_solution)
    # println("Temps CPU GRASP (s): ", grasp_cpu_time)

    return (heuristic_solution = solution_finale,
            heuristic_value = valeur_heur,
            cpu_time = cpu_time,
            glout_value = valeur_glout
            )
end


function utility(ensemble, C, A)
    return C[ensemble] / sum(A[:, ensemble])
end

function construction_gloutonne(C, A)
    n = length(C)
    solution_final = [] 
    solution = [] 
    # weight = 0
    elements = zeros(Bool, size(A, 1)) # suivi des éléments couverts
    # println(elements)

    # on calcule l'utilité de chaque ensemble
    for i in 1:n
        utility = C[i] / sum(A[:, i])
        push!(solution, (i, utility))
        # print("Ensemble $i, utilité : $utility\n")
    end
    # on trie par utilité décroissante
    sort!(solution, by = x -> x[2], rev = true)

    # on parcourt la liste des ensembles triés par utilité 
    while !isempty(solution)
        # on prend l'ensemble avec la meilleure utilité
        index, _ = popfirst!(solution)
        # on vérifie s'il y a conflit avec les éléments déjà couverts
        a_conflit = false
        for j in eachindex(elements)
            if A[j, index] == 1 && elements[j]
                a_conflit = true # il y a conflit avec un élément déjà couvert
                break
            end
        end

        # on l'ajoute à la solution seulement s'il n'y a pas de conflit
        if !a_conflit
            push!(solution_final, index)
            # on marque toutes ses contraintes comme prises
            for j in eachindex(elements)
                if A[j, index] == 1
                    elements[j] = true # on marque l'élément comme couvert
                    # weight += C[index] / 2
                end
            end
            # println("Ajout de l'ensemble $index, Elements couverts : $elements")
            # println("Poids actuel : $weight")
        end
    end

    return solution_final
end


function greedy_randomized_construction(C, A, α, rcl_size=10)
    n = length(C)
    m = size(A, 1)
    elements = zeros(Bool, m)  
    solution_final = []
    utilitys = []  

    # Initialiser les candidats
    for i in 1:n
        push!(utilitys, (i, utility(i, C, A)))
    end

    while !isempty(utilitys)
        u_min = minimum(x[2] for x in utilitys)
        u_max = maximum(x[2] for x in utilitys)
        uLimit = u_min + α * (u_max - u_min)
        # println("u_min: $u_min, u_max: $u_max, uLimit: $uLimit")

        rcl = filter(x -> x[2] >= uLimit, utilitys)
        # Trier par utilité décroissante
        sort!(rcl, by = x -> x[2], rev = true)

        # Réduire la RCL en prenant au plus rcl_size éléments (si rcl_size > length(rcl) on prend tout)
        if isempty(rcl)
            break
        end
        take_k = min(rcl_size, length(rcl))
        rcl = rcl[1:take_k]

        selected = rand(rcl)
        e = selected[1]

        push!(solution_final, e)

        # Marquer les contraintes couvertes
        for j in 1:m
            if A[j, e] == 1
                elements[j] = true
            end
        end

        # retirer ceux en conflit avec e
        new_utilitys = []
        for (cand, util) in utilitys
            if cand == e
                continue
            end

            # Vérifier si en conflit
            conflit = false
            for j in 1:m
                if A[j, cand] == 1 && elements[j]
                    conflit = true
                    break
                end
            end
            if !conflit
                push!(new_utilitys, (cand, util))
            end
        end
        utilitys = new_utilitys
    end

    return solution_final
end



function grasp(A,C, iteration)
    global_best = []
    global_best_value = -Inf

    for i in 1:iteration
        s = greedy_randomized_construction(C, A, 0.3)
        s = descente_locale(s, C, A)
        val = sum(C[j] for j in s)
        println("Itération $i, solution = ", s, "  valeur = ", val)
        if val > global_best_value
            global_best = copy(s)
            global_best_value = val
        end
    end
    println("GRASP meilleure valeur trouvée = ", global_best_value)
    return global_best
end


# Caches globaux pour optimiser les performances
const voisin_cache = Dict{Tuple{Int,Vararg{Int}}, Vector{Vector{Int}}}()
const evaluation_cache = Dict{Tuple{Int,Vararg{Int}}, Float64}()
const col_bitvectors_cache = Dict{UInt64, Tuple{Vector{BitVector}, Vector{Vector{Int}}}}()

# Fonction utilitaire pour évaluer une solution avec cache
function evaluer_solution(solution, C)
    if isempty(solution)
        return 0.0
    end
    key = Tuple(sort(solution))
    if haskey(evaluation_cache, key)
        return evaluation_cache[key]
    end
    val = sum(C[i] for i in solution)
    evaluation_cache[key] = val
    return val
end

# Pré-calcul optimisé des structures de données pour A
function get_col_structures(A)
    h = hash(A)
    if haskey(col_bitvectors_cache, h)
        return col_bitvectors_cache[h]
    end
    
    n = size(A, 2)
    colBits = Vector{BitVector}(undef, n)
    rows_by_col = Vector{Vector{Int}}(undef, n)
    
    for j in 1:n
        colBits[j] = BitVector(A[:, j] .== 1)
        rows_by_col[j] = findall(colBits[j])
    end
    
    result = (colBits, rows_by_col)
    col_bitvectors_cache[h] = result
    return result
end

function generer_voisinage(solution, C, A)
    # Vérifier le cache
    key = Tuple(sort(solution))
    if haskey(voisin_cache, key)
        return voisin_cache[key]
    end
    
    n = size(A, 2)
    m = size(A, 1)
    
    # Obtenir les structures pré-calculées
    colBits, rows_by_col = get_col_structures(A)
    
    # Set pour dédupliquer et ensemble de la solution courante
    voisins_set = Set{Vector{Int}}()
    solution_set = Set(solution)
    
    # Calculer la couverture actuelle avec cover_count
    cover_count = zeros(Int, m)
    for s in solution
        for r in rows_by_col[s]
            cover_count[r] += 1
        end
    end
    couverts = BitVector(cover_count .> 0)
    
    compteur11 = 0
    compteur21 = 0
    compteur01 = 0
    compteur10 = 0
    
    # 1-0 : Supprimer chaque élément (toujours valide)
    for s in solution
        new_sol = sort([x for x in solution if x != s])
        push!(voisins_set, new_sol)
        compteur10 += 1
    end
    
    # 0-1 : Ajouter un candidat si pas de conflit
    for cand in 1:n
        if cand in solution_set
            continue
        end
        if any(colBits[cand] .& couverts)
            continue
        end
        new_sol = sort(vcat(solution, cand))
        push!(voisins_set, new_sol)
        compteur01 += 1
    end
    
    # 1-1 : Échanger un élément avec un candidat
    for s in solution
        # Décrémenter temporairement la couverture
        for r in rows_by_col[s]
            cover_count[r] -= 1
        end
        couverts_temp = BitVector(cover_count .> 0)
        
        for cand in 1:n
            if cand in solution_set
                continue
            end
            if any(colBits[cand] .& couverts_temp)
                continue
            end
            new_sol = sort(vcat([x for x in solution if x != s], cand))
            push!(voisins_set, new_sol)
            compteur11 += 1
        end
        
        # Restaurer la couverture
        for r in rows_by_col[s]
            cover_count[r] += 1
        end
    end
    
    # 2-1 : Retirer 2 éléments et ajouter 1 candidat
    solution_arr = collect(solution)
    for i in 1:length(solution_arr)
        for j in (i+1):length(solution_arr)
            s1, s2 = solution_arr[i], solution_arr[j]
            
            # Décrémenter pour les deux
            for r in rows_by_col[s1]
                cover_count[r] -= 1
            end
            for r in rows_by_col[s2]
                cover_count[r] -= 1
            end
            couverts_temp = BitVector(cover_count .> 0)
            
            for cand in 1:n
                if cand in solution_set
                    continue
                end
                if any(colBits[cand] .& couverts_temp)
                    continue
                end
                new_sol = sort(vcat([x for x in solution if x != s1 && x != s2], cand))
                push!(voisins_set, new_sol)
                compteur21 += 1
            end
            
            # Restaurer
            for r in rows_by_col[s1]
                cover_count[r] += 1
            end
            for r in rows_by_col[s2]
                cover_count[r] += 1
            end
        end
    end
    
    voisins = collect(voisins_set)
    println("Généré: 0→1=$compteur01, 1→0=$compteur10, 1→1=$compteur11, 2→1=$compteur21, Total=$(length(voisins))")
    
    # Stocker dans le cache
    voisin_cache[key] = voisins
    return voisins
end


function descente_locale(solution_initiale, C, A)
    solution_courante = solution_initiale
    amelioration = true
    iteration = 0
    max_iterations = 1000  # Sécurité pour éviter boucles infinies
    
    # Ensemble des solutions visitées pour éviter les cycles
    visited = Set{Tuple{Int, Vararg{Int}}}()
    push!(visited, Tuple(sort(solution_courante)))

    while amelioration && iteration < max_iterations
        iteration += 1
        amelioration = false
        
        # Utiliser le cache pour évaluer la solution courante
        cout_courant = evaluer_solution(solution_courante, C)
        
        # Générer les voisins (utilise le cache automatiquement)
        voisins = generer_voisinage(solution_courante, C, A)
        println("Itération $iteration: $(length(voisins)) voisins, valeur actuelle = $cout_courant")
        
        # Trouver le meilleur voisin non visité
        meilleur_voisin = nothing
        meilleur_cout = cout_courant
        
        for voisin in voisins
            key_voisin = Tuple(sort(voisin))
            
            # Ignorer les solutions déjà visitées
            if key_voisin in visited
                continue
            end
            
            # Évaluer avec cache
            cout_voisin = evaluer_solution(voisin, C)
            
            if cout_voisin > meilleur_cout
                meilleur_voisin = voisin
                meilleur_cout = cout_voisin
            end
        end
        
        # Si on a trouvé une amélioration
        if meilleur_voisin !== nothing
            solution_courante = meilleur_voisin
            push!(visited, Tuple(sort(solution_courante)))
            amelioration = true
            println("  → Amélioration trouvée: nouvelle valeur = $meilleur_cout")
        else
            println("  → Optimum local atteint")
        end
    end
    
    if iteration >= max_iterations
        println("Limite d'itérations atteinte ($max_iterations)")
    end
    
    return solution_courante
end


# function experimentationSPP()
#     dossier = "Data"
#     fichiers = filter(f -> endswith(f, ".dat"), readdir(dossier, join=true))

#     resultats = DataFrame(
#         instance = String[],
#         my_solution_glouton = String[],
#         valeur_glouton = Float64[],
#         my_solution = String[],
#         my_z = Float64[],
#         my_cpu_time = Float64[]
#     )

#     for fichier in fichiers[1:min(10, length(fichiers))]  # max 10 instances
#         println("\n--- Traitement de l’instance : $(basename(fichier)) ---")

#         resultat = resoudreSPP(fichier)

#         # Récupération des données retournées par resoudreSPP()
#         solution_glouton = construction_gloutonne(loadSPP(fichier)...) 
#         valeur_glouton = resultat.glout_value
#         solution_finale = resultat.heuristic_solution
#         valeur_finale = resultat.heuristic_value
#         cpu_time = resultat.cpu_time

#         push!(resultats, (
#             instance = basename(fichier),
#             my_solution_glouton = join(solution_glouton, ","),
#             valeur_glouton = valeur_glouton,
#             my_solution = join(solution_finale, ","),
#             my_z = valeur_finale,
#             my_cpu_time = cpu_time
#         ))
#     end

#     # Sauvegarde CSV
#     CSV.write("resultats_SPP.csv", resultats)
# end


solution_heuristique = resoudreSPP(fname)