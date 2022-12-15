function compute_kkt_optimality(g::Vector{Float64}, on_bound::Vector{Int})
    length(g) == 0 && return 0.
    @views begin
        g_kkt = g .* on_bound
        free_set = (on_bound .== 0)
        g_kkt[free_set] = abs.(g[free_set])
    end
    return maximum(g_kkt)
end

function solve_heuristic_bvls(
    node::Node,
    A::Matrix{Float64},
    y::Vector{Float64},
    λ::Float64,
    M::Float64,
    bnbparams::BnbParams,
    )

    AS1 = A[:,node.branch.S1]
    nS1 = size(AS1)[2]
    lb = fill(-M,(nS1,))
    ub = fill(M,(nS1,))

    xS1 = zeros(Float64,nS1)
    on_bound = zeros(Int,nS1)
    free_set = (on_bound .== 0)
    active_set = .!free_set
    free_set = findall(free_set)

    r = AS1 * xS1 - y
    cost = r' * r
    g = AS1' * r
    cost_change = Inf

    # ----- Initialization loop ----- #

    while !isempty(free_set)

        AS1_free = AS1[:, free_set]
        y_free = y - AS1 * (xS1 .* active_set)
        z = AS1_free \ y_free
        lbv = (z .< lb[free_set])
        ubv = (z .> ub[free_set])
        v = lbv .| ubv

        if any(lbv)
            ind = free_set[lbv]
            xS1[ind] = lb[ind]
            active_set[ind] .= true
            on_bound[ind] .= -1
        end
        if any(ubv)
            ind = free_set[ubv]
            xS1[ind] = ub[ind]
            active_set[ind] .= true
            on_bound[ind] .= 1
        end

        ind = free_set[.!v]
        xS1[ind] = z[.!v]
        r = AS1 * xS1 - y
        g = AS1' * r
        any(v) || break
        free_set = free_set[.!v]
    end

    # ----- Main loop ----- #

    optimality = compute_kkt_optimality(g,on_bound)

    while true

        (optimality < bnbparams.tol) && break
        move_to_free = argmax(g .* on_bound)
        on_bound[move_to_free] = 0
        
        while true

            free_set = (on_bound .== 0)
            active_set = .!free_set
            free_set = findall(free_set)
            x_free = xS1[free_set]
            lb_free = lb[free_set]
            ub_free = ub[free_set]
            AS1_free = AS1[:,free_set]
            y_free = y - AS1 * (xS1 .* active_set)
            z = AS1_free \ y_free
            lbv = findall(z .< lb_free)
            ubv = findall(z .> ub_free)
            v = vcat(lbv,ubv)

            if !isempty(v)
                
                alphas = vcat(
                    lb_free[lbv] - x_free[lbv], 
                    ub_free[ubv] - x_free[ubv]
                ) ./ (z[v] - x_free[v])
                i = argmin(alphas)
                i_free = v[i]
                alpha  = alphas[i]
                x_free *= 1. - alpha
                x_free += alpha * z
                xS1[free_set] = x_free

                if i < length(lbv)
                    on_bound[free_set[i_free]] = -1
                else
                    on_bound[free_set[i_free]] = 1
                end
            else
                x_free = z
                xS1[free_set] = x_free
                break
            end
        end

        r = AS1 * xS1 - y
        cost_new = 0.5 * (r' * r)
        cost_change = cost - cost_new

        (cost_change < bnbparams.tol * cost) && break

        cost = cost_new
        g = AS1' * r
        optimality = compute_kkt_optimality(g,on_bound)
    end

    x = zeros(size(A)[2])
    x[node.branch.S1] = xS1
    obj = 0.5 * norm(r,2)^2 + λ * nS1
    
    return obj, x
end
