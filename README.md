# metaSPPstu: solving the Set Packing Problem (student version)

## Les différentes fonctions utilitaires

### resoudreSPP(fname)
Cette fonction prend en entré un fichier .dat et utilise reactive GRASP pour résoundre l'instance du SPP

### experimentationSPP()
Cette fonction lance 10 instances et lance sur chacune des instances l'heuristique gloutonne et fais une descente locale. Elle affiche dans la console les résultats obtenus.

### etude_parametres_grasp()
Cette fonction lance 10 instances et utilise GRASP avec différents paramètres (nombre d'itérations fix ou temps limité) et affiche les résultats dans la console.
exemple d'appel

# README simple — usage pour le professeur

Voici un README minimal qui liste les appels disponibles dans `src/livrableEI1.jl`.

Exemple d'usage (copier/colle dans le REPL ou un script) :

```julia
# ----- Instance à résoudre -----
fname = "dat/pb_1000rnd0300.dat"


#------ Appel Livrable1 ------
resoudreSPP(fname)
# experimentationSPP()

#------ Appel Livrable2 ------
# resoudreGRASP(fname, 0.4, 200)
# resoudreREACTIVEGRASP(fname, 200, 20)
# etude_grasp(mode="iterations", total_iterations=200)
# etude_reactivegrasp(mode="iterations", total_iterations=200)

#------ Appel Livrable3 ------
# resoudreAG(fname, 200, 500, 0.4)
# etude_AG(repeats=3)
```

Notes :
- Par défaut l'instance utilisée dans les exemples est `pb_200rnd0700.dat`.
- La fonction appelée par défaut (exécution simple) est `resoudreSPP(fname)`.
- Les autres appels sont fournis en commentaire : décommentez la ligne souhaitée pour l'exécuter.

Fichier principal : `src/livrableEI1.jl` — ouvrez-le pour voir la liste complète des fonctions et leurs signatures.

Si vous voulez que je fournisse une version encore plus brève ou différente (par ex. uniquement les signatures), dites‑le et j'ajuste.