function screen!(
    tree::Tree,
    node::Node,
    asviews::ActivesetViews,
    u::Vector{Float64},
    Atu::Vector{Float64},
    B::Matrix{Float64},
    A::Matrix{Float64},
    y::Vector{Float64},
    λ::Float64,
    M::Float64,
    )

    @views begin
        pivotSbar = @. M * (abs(Atu[node.branch.Sbar]) - λ/M)
        pivotSbar0 = max.(pivotSbar,0)
        pivotSbar1 = max.(-pivotSbar,0)
        pivotS1 = @. M * (abs(Atu[node.branch.S1]) - λ/M)
        dobj = 0.5 * (y'y - norm(y-u,2)^2) - sum(pivotSbar0) - sum(pivotS1) 
    end

    Sbar_to_S0 = Vector{Int}()
    Sbar_to_S1 = Vector{Int}()
    for (iSbar,i) in enumerate(node.branch.Sbar)
        if pivotSbar1[iSbar] > tree.ub - dobj
            push!(Sbar_to_S0,i)
            tree.scr0 += 1
        end
        if pivotSbar0[iSbar] > tree.ub - dobj
            push!(Sbar_to_S1,i)
            tree.scr1 += 1
        end
    end

    if !isempty(intersect(Sbar_to_S0,Sbar_to_S1))
        node.lb = Inf
        node.x = rand(size(A)[2])
        return true
    end

    for i in Sbar_to_S0
        fixto!(node,i,0)
        swap_Sbar_to_S0!(node,asviews,A,i)
    end

    for i in Sbar_to_S1
        fixto!(node,i,1)
        swap_Sbar_to_S1!(node,asviews,A,B,i)
    end

    return false
end
