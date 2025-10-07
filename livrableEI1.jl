
using JuMP
using GLPK
using LinearAlgebra
include("loadSPP.jl")
include("setSPP.jl")
include("getfname.jl")


# le ! dans les fonction pour indiquer que la fonction modifie ses arguments

function resoudreSPP(fname)
    C, A = loadSPP(fname)
    @show C
    @show A

    # contrôle du temps CPU
    t_start = time()

    # heuristique de construction gloutonne
    solution_initiale = construction_gloutonne(C, A)
    println("Solution initiale (gloutonne) : ", solution_initiale)

    # amélioration par descente locale
    solution_finale = descente_locale(solution_initiale, C, A)

    # valeur ^z (somme des coûts des ensembles sélectionnés)
    valeur_z = sum(C[solution_finale])
    println("Valeur ^z : ", valeur_z)

    t_end = time()
    println("Temps CPU écoulé : ", t_end - t_start, " secondes")

    return solution_finale
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
        utility = C[i] / sum(A[:, i])
        push!(solution, (i, utility))
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
                    weight += C[index] / 2
                end
            end
            # println("Ajout de l'ensemble $index, Elements couverts : $elements")
            println("Poids actuel : $weight")
        end
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

function generer_voisinage_1_1(solution, C, A)
   
    voisins = []
    n = size(A, 2) 
    
    # Pour chaque ensemble candidat non utilisé
    for candidat = 1:n
        if !(candidat in solution)
            # Détecter les conflits avec ce candidat
            conflits = detecter_conflits(candidat, solution, A)

            # Si exactement 1 conflit 1-1 exchange possible
            if length(conflits) == 1
                conflit = conflits[1]
                nouvelle_solution = setdiff(solution, [conflit]) # Retirer l'ensemble en conflit
                push!(nouvelle_solution, candidat)
                push!(voisins, nouvelle_solution)
            end
            # Si aucun conflit  0-1 exchange (ajout simple)
            if isempty(conflits)
                nouvelle_solution = copy(solution)
                push!(nouvelle_solution, candidat)
                push!(voisins, nouvelle_solution)
            end
        end
    end
    
    # Aussi essayer les 1-0 exchanges (retirer un ensemble)
    for ensemble in solution
        nouvelle_solution = setdiff(solution, [ensemble])
        if !isempty(nouvelle_solution)  # Éviter la solution vide
            push!(voisins, nouvelle_solution)
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
        voisins = generer_voisinage_1_1(solution_courante, C, A)
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




