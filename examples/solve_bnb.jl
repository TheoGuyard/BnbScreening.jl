include("../src/BnbScreening.jl")

# ---------------------------------------------------------------------------- #
# Example : Solving a L0-penalized least-square problem using a BnB algorithm  #
# ---------------------------------------------------------------------------- #

# Generate problem data (gaussian setup or toeplitz setup)
k, m, n, snr = 5, 500, 300, 10*log10(10)
A, y, λ, M = data_toeplitz(k,m,n,snr)
#k, m, n, snr = 5, 500, 1000, 10*log10(10)
#A, y, λ, M = data_gaussian(k,m,n,snr)

# BnB parameters
bnbparams = BnbParams(
    lb_method=ACTIVESET,
    # Relaxation solving method
    #   - ACTIVESET : the relaxation is solved using an Active-Set algorithm
    #   - LBMIP : the relaxation is solved using CPLEX convex solver
    ub_method=BVLS,
    # Heuristic method
    #   - BVLS : the heuristic is solved using the Bounded-Variable least-square algorithm
    #   - LBMIP : the heuristic is solved using CPLEX convex solver
    maxtime=1000,           # Maximum solution time allowed
    tol=1e-8,               # Integer tolerance, ie. |x|<=tol <=> x=0
    maxiter=1000,           # Maximum number of iterations allowed in the Active-Set method
    screening=false,        # Toogle branch-screening rules
    verbosity=true,         # Toogle verbosity
    showevery=1,            # If verbosity=true, number of nodes between two displays
)

# Solve the problem using the BnB algorithm with appropriate parameters
result_bnb = solve_bnb(A,y,λ,M,bnbparams=bnbparams)
