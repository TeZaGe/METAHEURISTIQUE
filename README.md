# metaSPPstu: solving the Set Packing Problem (student version)

## Les différentes fonctions utilitaires

### resoudreSPP(fname)
Cette fonction prend en entré un fichier .dat et utilise reactive GRASP pour résoundre l'instance du SPP

### experimentationSPP()
Cette fonction lance 10 instances et lance sur chacune des instances l'heuristique gloutonne et fais une descente locale. Elle affiche dans la console les résultats obtenus.

### etude_parametres_grasp()
Cette fonction lance 10 instances et utilise GRASP avec différents paramètres (nombre d'itérations fix ou temps limité) et affiche les résultats dans la console.
exemple d'appel
```julia
etude_parametres_grasp(mode="time", limit_time=60.0)
etude_parametres_grasp(mode="iteration", nb_iteration=200)
```

### etude_parametres_reactivegrasp()
Cette fonction lance 10 instances et utilise Reactive GRASP avec différents paramètres (nombre d'itérations fix ou temps limité) et affiche les résultats dans la console.
exemple d'appel
```julia
etude_parametres_reactivegrasp(mode="time", limit_time=60.0)
etude_parametres_reactivegrasp(mode="iteration", nb_iteration=200)
```

## Reste a faire 
- Finir CR_EI2.tex la partie graphique
- Implémentation de l'algorithme génétique 
- Faire CR_EI3.tex