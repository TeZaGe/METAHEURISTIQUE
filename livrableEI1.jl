
using JuMP
using GLPK
using LinearAlgebra

using CSV
using DataFrames
include("loadSPP.jl")
include("setSPP.jl")
include("getfname.jl")


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
    println("Temps CPU total (s): ", cpu_time)

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
    weight = 0
    elements = zeros(Bool, size(A, 1)) # suivi des éléments couverts
    # println(elements)

    # on calcule l'utilité de chaque ensemble
    for i in 1:n
        push!(solution, (i, utility(i, C, A)))
    end
    # on trie les ensembles par utilité décroissante
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


function greedy_randomized_construction(C, A, α=0.3, rcl_size=10)
    n = length(C)
    m = size(A, 1)
    elements = zeros(Bool, m)  # suivi des contraintes couvertes
    solution_final = []
    utilitys = []  # liste des (ensemble, utilité) disponibles

    # Initialiser les candidats
    for i in 1:n
        push!(utilitys, (i, utility(i, C, A)))
    end

    while !isempty(utilitys)
        
        if isempty(utilitys)
            break
        end
        u_min = minimum(x[2] for x in utilitys)
        u_max = maximum(x[2] for x in utilitys)
        uLimit = u_min + α * (u_max - u_min)

        
        rcl = filter(x -> x[2] >= uLimit, utilitys)

        
        # Trier par utilité décroissante et prendre les top rcl_size
        sort!(rcl, by = x -> x[2], rev = true)
        rcl = rcl[1:rcl_size]

        if isempty(rcl)
            break
        end

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


function detecter_conflits(ensemble_candidat, solution_actuelle, A)
    conflits = []
    
    # Pour chaque ensemble dans la solution actuelle
    for ensemble in solution_actuelle
        for j in 1:size(A, 1)
            if A[j, ensemble] == 1 && A[j, ensemble_candidat] == 1
                # Il y a un conflit alors on l'ajoute à la liste
                push!(conflits, ensemble)
                break
            end
        end
    end

    return conflits
end

function generer_voisinage(solution, C, A)
    voisins = Vector{Vector{Int}}()
    n = size(A, 2)
    m = size(A, 1)

    # Représentation booléenne des ensembles utilisés pour accès rapide
    utilisés = falses(n)
    for i in solution
        utilisés[i] = true
    end

    # Pré-calcul des lignes couvertes par la solution actuelle
    couverts = falses(m)
    for s in solution
        couverts .|= A[:, s] .== 1
    end

    # ---- Cas 1 : Ajout (0→1) ----
    for cand in 1:n
        if !utilisés[cand]
            # Vérifie qu’il n’y a pas de conflit
            conflit = any((A[:, cand] .== 1) .& couverts)
            if !conflit
                push!(voisins, vcat(solution, [cand]))
            end
        end
    end

    # ---- Cas 2 : Retrait (1→0) ----
    for s in solution
        new_sol = setdiff(solution, [s])
        push!(voisins, new_sol)
    end

    # ---- Cas 3 : Échange (1→1) ----
    for s in solution
        # On retire temporairement s
        couverts_temp = copy(couverts)
        couverts_temp .&= .!(A[:, s] .== 1)

        for cand in 1:n
            if !utilisés[cand]
                # on vérifie que cand ne crée pas de conflit avec le reste
                conflit = any((A[:, cand] .== 1) .& couverts_temp)
                if !conflit
                    new_sol = setdiff(solution, [s])
                    push!(new_sol, cand)
                    push!(voisins, new_sol)
                end
            end
        end
    end

    return voisins
end
function descente_locale(solution_initiale, C, A)
    # on prend la solution initiale
    solution_courante = solution_initiale
    amelioration = true

    # on continue tant qu'on trouve une amélioration
    while amelioration
        amelioration = false
        voisins = generer_voisinage(solution_courante, C, A)
        # on évalue chaque voisin
        for voisin in voisins
            # on calcule le coût de la solution courante et du voisin
            cout_courant = sum(C[solution_courante])
            cout_voisin = sum(C[voisin])

            if cout_voisin > cout_courant
                solution_courante = voisin
                amelioration = true
                break
            end
        end
    end

    return solution_courante
end

function experimentationSPP()
    dossier = "Data"
    fichiers = filter(f -> endswith(f, ".dat"), readdir(dossier, join=true))

    resultats = DataFrame(
        instance = String[],
        my_solution_glouton = String[],
        valeur_glouton = Float64[],
        my_solution = String[],
        my_z = Float64[],
        my_cpu_time = Float64[]
    )

    for fichier in fichiers[1:min(10, length(fichiers))]  # max 10 instances
        println("\n--- Traitement de l’instance : $(basename(fichier)) ---")

        resultat = resoudreSPP(fichier)

        # Récupération des données retournées par resoudreSPP()
        solution_glouton = construction_gloutonne(loadSPP(fichier)...) 
        valeur_glouton = resultat.glout_value
        solution_finale = resultat.heuristic_solution
        valeur_finale = resultat.heuristic_value
        cpu_time = resultat.cpu_time

        push!(resultats, (
            instance = basename(fichier),
            my_solution_glouton = join(solution_glouton, ","),
            valeur_glouton = valeur_glouton,
            my_solution = join(solution_finale, ","),
            my_z = valeur_finale,
            my_cpu_time = cpu_time
        ))
    end

    # Sauvegarde CSV
    CSV.write("resultats_SPP.csv", resultats)
end

