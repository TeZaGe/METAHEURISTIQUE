# =========================================================================== #
# Compliant julia 1.x

# Using the following packages
using JuMP, GLPK
using LinearAlgebra

include("loadSPP.jl")
include("setSPP.jl")
include("getfname.jl")
include("livrableEI1.jl")

# =========================================================================== #

# Loading a SPP instance
println("\nLoading...")
fname = "Data/pb_100rnd0100.dat"
C, A = loadSPP(fname)
# @show C
# @show A


# Solving a SPP instance with GLPK
# println("\nSolving...")
# t_start = time()
# solverSelected = GLPK.Optimizer
# spp = setSPP(C, A)

# set_optimizer(spp, solverSelected)
# optimize!(spp)

# t_end = time()
# cpu_time = t_end - t_start
# println("Temps CPUt optimizer (s): ", cpu_time)

# Displaying the results
# println("z = ", objective_value(spp))
# print("x = "); println(value.(spp[:x]))
# =========================================================================== #
# solution_heuristique = resoudreSPP(fname)

experimentationSPP()
