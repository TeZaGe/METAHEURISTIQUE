# 📦 Métaheuristiques : Résolution du problème SPP

Ce projet regroupe plusieurs métaheuristique (Grasp, Algo génétique ...) pour le problème de Set Packing (SPP). Le point d'entrée est le fichier src/LivrableEI1.jl.

## 🚀 Guide d'utilisation

1. Choisir l'instance

Modifiez la variable fname a la fin du fichier src/livrableEI1.jl pour changer le jeu de données :
```julia
fname = "dat/pb_200rnd0700.dat"
```


2. Sélectionner l'algorithme

Le fichier src/livrableEI1.jl est structuré par "Livrables". Par défaut, seule l'heuristique gloutonne (resoudreSPP) est active.

Pour utiliser les méthodes avancées (GRASP, Génétique), il suffit de décommenter les lignes correspondantes (enlever le #) et de commenter celles que vous ne voulez pas exécuter.

Pour les différentes fonction d'étude statistique, les instances sont déjà sélectionnées dans les fonctions.

## 🟢 Livrable 1 : Heuristique de construction (Activé par défaut)

```julia
#------ Appel Livrable 1 ------
resoudreSPP(fname)        # Résolution gloutonne simple
# experimentationSPP()    # Tests statistiques 
```

## 🟠 Livrable 2 : GRASP & Reactive GRASP

Pour activer, décommentez les lignes ci-dessous :

les fonctions etude_grasp  et e possède deux modes : "iterations" ou on spécifie un nombre total d'itérations, et "time" ou on spécifie un temps total d'exécution en secondes.

```julia
#------ Appel Livrable 2 ------
# resoudreGRASP(fname, 0.4, 200)           # GRASP (alpha, itérations)
# resoudreREACRIVEGRASP(fname, 200, 20)    # Reactive GRASP (itérations, bloc_size)
# etude_grasp(mode="iterations", total_iterations=200) # Stats GRASP 
# etude_reactive_grasp(mode="iterations", total_iterations=200) # Stats Reactive GRASP
```

## 🔵 Livrable 3 : Algorithme Génétique

Pour activer, décommentez les lignes ci-dessous :
Pour
```julia
#------ Appel Livrable 3 ------
# resoudreAG(fname, 200, 500, 0.4)  # Algo Génétique (pop_size, gen, cross_rate)
# etude_AG(repeats=3)               # Stats Génétique 
```

## 📂 Organisation des dossiers

src/ : Code source.

dat/ : Fichiers d'instances de données.

res/ : Résultats de sortie.

doc/ : Documentation et graphiques.

## Lancement

Depuis la racine du projet :
```shell
julia src/livrableEI1.jl
```
ou directement dans julia :
```julia
 include("src/livrableEI1.jl")
```


