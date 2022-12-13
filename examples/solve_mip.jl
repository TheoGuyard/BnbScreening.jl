include("../src/BnbScreening.jl")

# ---------------------------------------------------------------------------- #
# Example : Solving a L0-penalized least-square problem using a MIP solver     #
# ---------------------------------------------------------------------------- #

# Generate problem data (gaussian setup or toeplitz setup)
k, m, n, snr = 5, 500, 300, 10*log10(10)
A, y, λ, M = data_toeplitz(k,m,n,snr)
#k, m, n, snr = 5, 500, 1000, 10*log10(10)
#A, y, λ, M = data_gaussian(k,m,n,snr)

# MIP solver parameters
mipparams = MipParams(
    maxtime = 1000,         # Maximum solution time in seconds
    silent = false,         # Toogle solver verbosity
)

# Solve the problem using a MIP solver with appropriate parameters
result_mip = solve_mip(A,y,λ,M,mipparams=mipparams)

# Display solver results and statistics
println("\nMIP")
println("  Status    : $(result_mip.termination_status)")
println("  Objective : $(result_mip.objective_value)")
println("  Nodes     : $(result_mip.node_count)")
println("  Timer     : $(result_mip.solve_time)")
