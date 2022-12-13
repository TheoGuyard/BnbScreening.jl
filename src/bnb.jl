function best_lb(tree::Tree)
    isempty(tree.active) && return 0.
    return maximum([node.lb for node in tree.active])
end

function branch!(tree::Tree, node::Node)
    prune(tree,node) && return nothing
    j = branching_index(node)
    isa(j,Nothing) && return nothing
    node_j0 = Node(node,j,0)
    node_j1 = Node(node,j,1)
    push!(tree.queue,node_j0)
    push!(tree.queue,node_j1)
end

function branching_index(node::Node)
    isempty(node.branch.Sbar) && return nothing
    return node.branch.Sbar[argmax(abs.(node.x[node.branch.Sbar]))]
end

function depth(node::Node)
    depth = 0
    while !isa(node,Nothing)
        depth += 1
        node = node.parent
    end
    return depth
end

function elapsed_time(tree::Tree)
    return (Dates.now() - tree.start_time).value / 1000.
end

function fixto!(node::Node, j::Int, jval::Int)
    jpos = findfirst(i->i==j,node.branch.Sbar)
    isa(jpos,Nothing) && error("Branching index is already fixed")
    if jval == 0
        push!(node.branch.S0,popat!(node.branch.Sbar,jpos))
        node.x[j] = 0.
    elseif jval == 1
        push!(node.branch.S1,popat!(node.branch.Sbar,jpos))
    end
    return nothing
end

function global_gap(tree::Tree)
    return (tree.ub - best_lb(tree)) / tree.ub
end

function has_perfect_relaxation(node::Node, M::Float64, bnbparams::BnbParams)
    return !(any(bnbparams.tol .< abs.(node.x[node.branch.Sbar]) .< M - bnbparams.tol))
end

function lower_bound!(
    tree::Tree,
    node::Node,
    A::Matrix{Float64},
    y::Vector{Float64},
    λ::Float64,
    M::Float64,
    bnbparams::BnbParams
    )
    if bnbparams.lb_method == ACTIVESET
        solve_relax_activeset!(tree,node,A,y,λ,M,bnbparams)
    elseif bnbparams.lb_method == LBMIP
        solve_relax_mip!(node,A,y,λ,M)
    else
        error("Paramter `bnbparams.lb_method` not understood.")
    end
    return nothing
end

function next_node!(tree::Tree)
    return pop!(tree.queue)
end

function processed!(tree::Tree, node::Node)
    if !prune(tree,node)
        push!(tree.active,node)
    end
    tree.nexpl += 1
    return nothing
end

function prune(tree::Tree, node::Node)
    return (tree.ub < node.lb)
end

function terminated!(tree::Tree, bnbparams::BnbParams)
    if isempty(tree.queue)
        tree.termination_status = MOI.OPTIMAL
    elseif elapsed_time(tree) > bnbparams.maxtime
        tree.termination_status = MOI.TIME_LIMIT
    end
    return (tree.termination_status != MOI.OPTIMIZE_NOT_CALLED)
end

function update_upper_bound!(tree::Tree, ub::Float64, inc::Vector{Float64})
    if ub < tree.ub
        tree.ub = ub
        tree.inc = inc
        return true
    end
    return false
end

function update_active_nodes!(tree::Tree)
    filter!(node->!prune(tree,node),tree.active)
    return nothing
end

function upper_bound!(
    tree::Tree,
    node::Node,
    A::Matrix{Float64},
    y::Vector{Float64},
    λ::Float64,
    M::Float64,
    bnbparams::BnbParams
    )

    perfect = false

    if isempty(node.branch.S1)
        ub = 0.5*y'*y
        inc = zeros(Float64,size(A)[2])
    elseif has_perfect_relaxation(node,M,bnbparams)
        ub = node.lb
        inc = node.x
        perfect = true
    else
        if bnbparams.ub_method == BVLS
            ub, inc = solve_heuristic_bvls(node,A,y,λ,M,bnbparams)
        elseif bnbparams.ub_method == UBMIP
            ub, inc = solve_heuristic_mip(node,A,y,λ,M)
        else
            error("Paramter `bnbparams.ub_method` not understood.")
        end
    end
    updated = update_upper_bound!(tree,ub,inc)
    updated && update_active_nodes!(tree)
    return perfect
end

# ============================================================================ #

"""
    solve_bnb(
        A::Matrix{Float64}, 
        y::Vector{Float64}, 
        λ::Float64, 
        M::Float64;
        bnbparams::BnbParams=BnbParams(),
        )

Solve the L0-penalized least-square problem 

    min_x   (1/2) norm(y-Ax,2)^2 + λ norm(x,0)  
    st.     norm(x,Inf) <= M

using a Branch-and-Bound algorithm. This algorithm can use branch-screening 
tests to enhance its efficiency. Returns a [`@BnbResults`](@ref) struct 
containing the optimization results and statistics. See [`@BnbParams`](@ref) for 
informations about additional method parameters.
"""
function solve_bnb(
    A::Matrix{Float64},
    y::Vector{Float64},
    λ::Float64,
    M::Float64;
    bnbparams::BnbParams=BnbParams(),
    )

    validate_data(A,y,λ,M)
    tree = Tree(A,y)
    bnbparams.verbosity && display_head()

    while !terminated!(tree,bnbparams)
        node = next_node!(tree)
        lower_bound!(tree,node,A,y,λ,M,bnbparams)
        perfect = upper_bound!(tree,node,A,y,λ,M,bnbparams)
        processed!(tree,node)
        perfect || branch!(tree,node)
        bnbparams.verbosity && display_node(tree,node)
    end

    results = BnbResults(tree)
    bnbparams.verbosity && display_tail(results)
    
    return results
end
