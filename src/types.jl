# ========================================== #
# Active-Set support configuration and views #
# ========================================== #

mutable struct Support
    Sbar0::Vector{Int}      # Indices verifying 0 = xi      with i in Sbar
    Sbarin::Vector{Int}     # Indices verifying 0 < xi < M  with i in Sbar
    Sbarbox::Vector{Int}    # Indices verifying     xi = M  with i in Sbar
    S1in::Vector{Int}       # Indices verifying     xi < M  with i in S1
    S1box::Vector{Int}      # Indices verifying     xi = M  with i in S1
end

mutable struct ActivesetViews
    # Current iterate
    xSbarin::Vector{Float64}
    xSbarbox::Vector{Float64}
    xS1in::Vector{Float64}
    xS1box::Vector{Float64}

    # Sign of x
    sSbarin::Vector{Float64}
    sSbarbox::Vector{Float64}
    sS1in::Vector{Float64}
    sS1box::Vector{Float64}

    # New iterate
    xnewSbarin::Vector{Float64}
    xnewS1in::Vector{Float64}

    # Submatrices
    ASbarin::Matrix{Float64}
    ASbarbox::Matrix{Float64}
    AS1in::Matrix{Float64}
    AS1box::Matrix{Float64}
end


# =========================== #
# Branch-and-Bound components #
# =========================== #

mutable struct Branch
    S0::Vector{Int}     # Indices fixed to zero
    S1::Vector{Int}     # Indices fixed to non-zero
    Sbar::Vector{Int}   # Indices not fixed
end

mutable struct Node
    parent::Union{Nothing,Node}     # Parent node
    branch::Branch                  # Current branch (fixed variables)
    branch_on::Tuple{Int,Int}       # Branching index above node
    x::Vector{Float64}              # Relaxation solution
    lb::Float64                     # Relaxation objective value

    # Active set data cashed in each node
    support::Support
    GS1in_inv::Matrix{Float64}      
    GSbarin_inv::Matrix{Float64}
end

mutable struct Tree
    termination_status::MOI.TerminationStatusCode   # Optimization status
    start_time::DateTime                            # Start time
    nexpl::Int                                      # Number of nodes explored
    queue::Vector{Node}                             # Nodes to explore
    active::Vector{Node}                            # Active nodes (with node.lb < tree.ub)
    ub::Float64                                     # Best upper bound
    inc::Vector{Float64}                            # Incumbent solution
    scr0::Int                                       # Indices fixed to 0 with screening tests
    scr1::Int                                       # Indices fixed to 1 with screening tests
end


# ======================================= #
# Branch-and-Bound parameters and outputs #
# ======================================= #

# Available solving methods for the relaxed problem (lower bounding step)
@enum LbMethod begin
    LBMIP       # MIP solver
    ACTIVESET   # Active-Set algorihtm
end

# Available heuristic methods (upper bounding step)
@enum UbMethod begin
    UBMIP   # MIP solver
    BVLS    # Bounded Variable Least-Square algorithm
end

struct BnbParams
    lb_method::LbMethod     # Method used to solve node relaxations
    ub_method::UbMethod     # Heuristic to compute node uppper bounds
    maxtime::Int            # Maximum solving time in seconds
    tol::Float64            # Integer tolerance
    maxiter::Int            # Max. iter. of the relaxation solution algorithm
    screening::Bool         # Toogle screening
    verbosity::Bool         # Whether to enable iterations displays
    showevery::Int          # Node display interval
end

struct BnbResults
    termination_status::MOI.TerminationStatusCode   # Optimization status
    solve_time::Float64                             # Solution time in seconds
    node_count::Int                                 # Number of nodes explored
    objective_value::Float64                        # Objective value (if optimal)
    x::Vector{Float64}                              # Optimizer (if optimal)
    scr0::Int                                       # Number of screning tests passed (xi = 0)
    scr1::Int                                       # Number of screning tests passed (xi != 0)
end


# ================================= #
# MIP solver parameters and outputs #
# ================================= #

struct MipParams
    maxtime::Int    # Maximum solution time in seconds
    silent::Bool    # Set solver silent
end

struct MipResults
    termination_status::MOI.TerminationStatusCode   # Optimization status
    solve_time::Float64                             # Solution time in seconds
    node_count::Int                                 # Number of nodes explored
    objective_value::Float64                        # Objective value (if optimal)
    x::Vector{Float64}                              # Optimizer (if optimal)
end


# ============ #
# Constructors #
# ============ #

function Support(n::Int)
    return Support(
        collect(1:n), 
        Vector{Int}(), 
        Vector{Int}(), 
        Vector{Int}(), 
        Vector{Int}()
    )
end

function ActivesetViews(node::Node, A::Matrix{Float64}, bnbparams::BnbParams)
    xSbarin = node.x[node.support.Sbarin]
    xSbarbox = node.x[node.support.Sbarbox]
    xS1in = node.x[node.support.S1in]
    xS1box = node.x[node.support.S1box]
    sSbarin = tolsign.(xSbarin, bnbparams.tol)
    sSbarbox = tolsign.(xSbarbox, bnbparams.tol)
    sS1in = tolsign.(xS1in, bnbparams.tol)
    sS1box = tolsign.(xS1box, bnbparams.tol)
    xnewSbarin = copy(xSbarin)
    xnewS1in = copy(xS1in)
    ASbarin = A[:,node.support.Sbarin]
    ASbarbox = A[:,node.support.Sbarbox]
    AS1in = A[:,node.support.S1in]
    AS1box = A[:,node.support.S1box]
    return ActivesetViews(
        xSbarin,
        xSbarbox,
        xS1in,
        xS1box,
        sSbarin,
        sSbarbox,
        sS1in,
        sS1box,
        xnewSbarin,
        xnewS1in,
        ASbarin,
        ASbarbox,
        AS1in,
        AS1box,
    )
end

function Branch(n::Int)
    return Branch(Vector{Int}(), Vector{Int}(), collect(1:n))
end

function Node(A::Matrix{Float64})
    m, n = size(A)
    return Node(
        nothing, 
        Branch(n), 
        (-1,-1), 
        zeros(Float64,n),
        0.,
        Support(n),
        zeros(Float64,0,0),
        zeros(Float64,0,0),
        )
end

function Node(parent::Node, j::Int, jval::Int)
    node = Node(
        parent,
        copy(parent.branch),
        (j,jval),
        copy(parent.x),
        parent.lb,
        copy(parent.support),
        copy(parent.GS1in_inv),
        copy(parent.GSbarin_inv),
        )
    fixto!(node,j,jval)
    return node
end

function Tree(A::Array{Float64,2}, y::Vector{Float64})
    return Tree(
        MOI.OPTIMIZE_NOT_CALLED,
        Dates.now(),
        0,
        Vector{Node}([Node(A)]),
        Vector{Node}(),
        0.5 * norm(y, 2)^2,
        zeros(Float64,size(A)[2]),
        0,
        0,
    )
end

function BnbParams(;
    lb_method=ACTIVESET,
    ub_method=BVLS,
    maxtime=1000,
    tol=1e-8,
    maxiter=1000,
    screening=false,
    verbosity=true,
    showevery=1,
    )
    maxtime > 0 || error("Parameter `maxtime` must be positive")
    tol > 0 || error("Parameter `tol` must be positive")
    maxiter > 0 || error("Parameter `maxiter` must be positive")
    showevery > 0 || error("Parameter `showevery` must be positive")
    return BnbParams(
        lb_method,
        ub_method,
        maxtime,
        tol,
        maxiter,
        screening,
        verbosity,
        showevery
        )
end

function BnbResults(tree::Tree)
    return BnbResults(
        tree.termination_status,
        elapsed_time(tree),
        tree.nexpl,
        tree.ub,
        tree.inc,
        tree.scr0,
        tree.scr1,
    )
end

function MipParams(;maxtime=1000,silent=true)
    maxtime > 0 || error("Parameter `maxtime` must be positive")
    return MipParams(maxtime, silent)
end

function MipResults(problem::JuMP.Model)
    return MipResults(
        termination_status(problem),
        solve_time(problem),
        convert(Int,node_count(problem)),
        objective_value(problem),
        value.(problem[:x])
    )
end


# ====================== #
# Base module extensions #
# ====================== #


function Base.copy(support::Support)
    return Support(
        copy(support.Sbar0), 
        copy(support.Sbarin), 
        copy(support.Sbarbox), 
        copy(support.S1in), 
        copy(support.S1box)
    )
end

function Base.copy(branch::Branch)
    return Branch(copy(branch.S0), copy(branch.S1), copy(branch.Sbar))
end
