# ICASSP experient - Performances

include("../../src/BnbScreening.jl")
using StatsBase

################################ Configuration #################################

# This is the only part of code you may have to edit

setup = "toeplitz"          # Setup to use : "gaussian" or "toeplitz"
k = 5                       # Sparsity level
m, n = 500, 300            # Dictionnary dimension
snr = 10. * log10(10.)      # SNR of the observation noise
maxtime = 1000              # Maximum solution time (in second) allowed
maxiter = 1000              # Maximum Active-Set iterations allowed per relaxation resolution
tol = 1.e-8                 # Integer tolerance (ie, |x|<tol <=> x=0)
repeats = 10               # Number of repeats in the experiment


################################## Experiment ##################################

function print_and_log(s::String, f::IOStream)
    print(s)
    write(f,s)
end

function precompilation(
    setup::String,
    k::Int,
    m::Int,
    n::Int,
    snr::Float64,
    maxtime::Int,
    maxiter::Int,
    tol::Float64,
    )
    print_and_log("Precompiling methods ...\n",logfile)
    try
        if setup == "gaussian"
            A, y, λ, M = data_gaussian(k,m,n,snr)
        elseif setup == "toeplitz"
            A, y, λ, M = data_toeplitz(k,m,n,snr)
        else
            error("Parameter `setup` must be either \"gaussian\" or \"toeplitz\".")
        end
        result_mip = solve_mip(A,y,λ,M,mipparams=MipParams(maxtime=maxtime))
        result_bnb = solve_bnb(A,y,λ,M,bnbparams=BnbParams(maxtime=maxtime,maxiter=maxiter,tol=tol,verbosity=false))
        result_bnbs = solve_bnb(A,y,λ,M,bnbparams=BnbParams(maxtime=maxtime,maxiter=maxiter,tol=tol,verbosity=false,screening=true))
    catch e
        print_and_log("Precompilation failed with the following error\n",logfile)
        print_and_log("$e\n",logfile)
        print_and_log("Stopping experiment\n",logfile)
        rethrow(e)
    end
    print_and_log("\n",logfile)
    return nothing
end

logfile_name = string(
    @__DIR__,
    "/",
    Dates.format(Dates.now(),"yyyy-mm-dd-HH:MM:SS"),
    ".txt"
    )
logfile = open(logfile_name, "w")

print_and_log("--------------------------------\n",logfile)
print_and_log("ICASSP experiment - Performances\n\n",logfile)
print_and_log("Configuration\n",logfile)
print_and_log("  setup   : $(setup)\n",logfile)
print_and_log("  k       : $(k)\n",logfile)
print_and_log("  m, n    : $(m), $(n)\n",logfile)
print_and_log("  snr     : $(snr)\n",logfile)
print_and_log("  maxtime : $(maxtime)\n",logfile)
print_and_log("  maxiter : $(maxiter)\n",logfile)
print_and_log("  tol     : $(tol)\n",logfile)
print_and_log("  repeats : $(repeats)\n",logfile)
print_and_log("\n",logfile)

fail_mip = Vector{Union{Bool,Missing}}(missing,repeats)
fail_bnb = Vector{Union{Bool,Missing}}(missing,repeats)
fail_bnbs = Vector{Union{Bool,Missing}}(missing,repeats)
node_mip = Vector{Union{Int64,Missing}}(missing,repeats)
node_bnb = Vector{Union{Int64,Missing}}(missing,repeats)
node_bnbs = Vector{Union{Int64,Missing}}(missing,repeats)
time_mip = Vector{Union{Float64,Missing}}(missing,repeats)
time_bnb = Vector{Union{Float64,Missing}}(missing,repeats)
time_bnbs = Vector{Union{Float64,Missing}}(missing,repeats)

precompilation(setup,k,m,n,snr,maxtime,maxiter,tol)

for i in 1:repeats
    print_and_log("Repeat $i/$repeats\n",logfile)
    try
        # Generate data
        if setup == "gaussian"
            A, y, λ, M = data_gaussian(k,m,n,snr)
        elseif setup == "toeplitz"
            A, y, λ, M = data_toeplitz(k,m,n,snr)
        end

        # Run the tree methods
        result_mip = solve_mip(A,y,λ,M,mipparams=MipParams(maxtime=maxtime))
        result_bnb = solve_bnb(A,y,λ,M,bnbparams=BnbParams(maxtime=maxtime,maxiter=maxiter,tol=tol,verbosity=false))
        result_bnbs = solve_bnb(A,y,λ,M,bnbparams=BnbParams(maxtime=maxtime,maxiter=maxiter,tol=tol,verbosity=false,screening=true))

        # Gather results
        if result_mip.termination_status == MOI.OPTIMAL
            fail_mip[i] = false
            node_mip[i] = result_mip.node_count
            time_mip[i] = result_mip.solve_time
        else
            fail_mip[i] = true
        end
        if result_bnb.termination_status == MOI.OPTIMAL
            fail_bnb[i] = false
            node_bnb[i] = result_bnb.node_count
            time_bnb[i] = result_bnb.solve_time
        else
            fail_bnb[i] = true
        end
        if result_bnbs.termination_status == MOI.OPTIMAL
            fail_bnbs[i] = false
            node_bnbs[i] = result_bnbs.node_count
            time_bnbs[i] = result_bnbs.solve_time
        else
            fail_bnbs[i] = true
        end

        # Log results
        print_and_log("  Fails\n",logfile)
        print_and_log("    MIP  : $(fail_mip[i])\n",logfile)
        print_and_log("    BnB  : $(fail_bnb[i])\n",logfile)
        print_and_log("    SBnB : $(fail_bnbs[i])\n",logfile)
        print_and_log("  Node count\n",logfile)
        print_and_log("    MIP  : $(node_mip[i])\n",logfile)
        print_and_log("    BnB  : $(node_bnb[i])\n",logfile)
        print_and_log("    SBnB : $(node_bnbs[i])\n",logfile)
        print_and_log("  Solve time\n",logfile)
        print_and_log("    MIP  : $(round(time_mip[i],digits=3)) seconds\n",logfile)
        print_and_log("    BnB  : $(round(time_bnb[i],digits=3)) seconds\n",logfile)
        print_and_log("    SBnB : $(round(time_bnbs[i],digits=3)) seconds\n",logfile)

    catch e
        if e == InterruptException()
            close(logfile)
            rethrow(e)
        end
        print_and_log("warning : Repeat $i failed with the following error\n",logfile)
        print_and_log("$e\n",logfile)
        print_and_log("Skipping this repeat\n",logfile)
    end
end


################################### Results ####################################

print_and_log("\n",logfile)
print_and_log("Results\n",logfile)
print_and_log("  Fails (mean over $repeats repeats)\n",logfile)
print_and_log("    MIP  : $(100 * mean(skipmissing(fail_mip)))%\n",logfile)
print_and_log("    BnB  : $(100 * mean(skipmissing(fail_bnb)))%\n",logfile)
print_and_log("    SBnB : $(100 * mean(skipmissing(fail_bnbs)))%\n",logfile)
print_and_log("  Node count (mean over $repeats repeats)\n",logfile)
print_and_log("    MIP  : $(mean(skipmissing(node_mip)))\n",logfile)
print_and_log("    BnB  : $(mean(skipmissing(node_bnb)))\n",logfile)
print_and_log("    SBnB : $(mean(skipmissing(node_bnbs)))\n",logfile)
print_and_log("  Solve time (mean over $repeats repeats)\n",logfile)
print_and_log("    MIP  : $(round(mean(skipmissing(time_mip)),digits=3)) seconds\n",logfile)
print_and_log("    BnB  : $(round(mean(skipmissing(time_bnb)),digits=3)) seconds\n",logfile)
print_and_log("    SBnB : $(round(mean(skipmissing(time_bnbs)),digits=3)) seconds\n",logfile)
print_and_log("--------------------------------",logfile)
close(logfile)
