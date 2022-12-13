function solve_relax_mip!(
    node::Node,
    A::Matrix{Float64},
    y::Vector{Float64},
    λ::Float64,
    M::Float64,
    )
    
    relax = Model(CPLEX.Optimizer)
    set_silent(relax)

    gram = A'*A
    Aty = A'*y
    
    n = size(A)[2]
    @variable(relax, x[1:n])
    @variable(relax, xabs[1:n])
    @objective(relax, Min, 
        0.5*(y'*y) - Aty'*x + (0.5*x)'*(gram*x) + 
        (λ/M)*sum(xabs[node.branch.Sbar]) + 
        λ*length(node.branch.S1)
    )
    @constraint(relax, box, -M .<= x .<= M)
    @constraint(relax, abs_pos, xabs .>= x)
    @constraint(relax, abs_neg, xabs .>= -x)
    @constraint(relax, setzero, x[node.branch.S0] .== 0)
    @constraint(relax, abs_setzero, xabs[node.branch.S0] .== 0)

    for i in 1:n
        set_start_value(x[i], node.x[i])
        set_start_value(xabs[i], abs(node.x[i]))
    end

    optimize!(relax)

    node.lb = objective_value(relax)
    node.x = value.(x)
end

function solve_heuristic_mip(
    node::Node,
    A::Matrix{Float64},
    y::Vector{Float64},
    λ::Float64,
    M::Float64,
    )

    heur = Model(CPLEX.Optimizer)
    set_silent(heur)

    nS1 = length(node.branch.S1)
    AS1 = A[:,node.branch.S1]
    xS1 = node.x[node.branch.S1]
    
    @variable(heur, x[1:nS1])
    @objective(heur, Min, 0.5*(y-AS1*x)'*(y-AS1*x)+ λ*length(node.branch.S1))
    @constraint(heur, box, -M .<= x .<= M)

    for i in 1:nS1
        set_start_value(x[i], xS1[i])
    end

    optimize!(heur)

    inc = zeros(Float64,size(A)[2])
    inc[node.branch.S1] = value.(x)

    return objective_value(heur), inc
end

# ============================================================================ #

"""
    solve_mip(
        A::Matrix{Float64}, 
        y::Vector{Float64}, 
        λ::Float64, 
        M::Float64;
        mipparams::MipParams=MipParams(),
        )

Solve the L0-penalized least-square problem 

    min_x   (1/2) norm(y-Ax,2)^2 + λ norm(x,0)  
    st.     norm(x,Inf) <= M

using a MIP solver (CPLEX). Returns a [`@MipResults`](@ref) struct containing
the optimization results and statistics. See [`@MipParams`](@ref) for 
informations about additional method parameters.
"""
function solve_mip(
    A::Matrix{Float64}, 
    y::Vector{Float64}, 
    λ::Float64, 
    M::Float64;
    mipparams::MipParams=MipParams(),
    )

    problem = Model(CPLEX.Optimizer)
    set_optimizer_attribute(problem,"CPX_PARAM_TILIM",mipparams.maxtime)
    mipparams.silent && set_silent(problem)

    gram = A'*A
    Aty  = A'*y
    
    n = size(A)[2]
    @variable(problem, x[1:n])
    @variable(problem, z[1:n], Bin)
    @objective(problem, Min, 0.5*(y'*y) - Aty'*x + (0.5*x)'*(gram*x) + λ*sum(z))
    @constraint(problem, box_lower, -M.*z .<= x)
    @constraint(problem, box_upper, x .<= M.*z)

    optimize!(problem)
    result = MipResults(problem)

    return result
end
