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
    n = size(AS1)[2]
    lb = fill(-M,(n,))
    ub = fill(M,(n,))

    x = zeros(Float64,n)
    on_bound = zeros(Int,n)
    free_set = (on_bound .== 0)
    active_set = .!free_set
    free_set = findall(free_set)

    r = AS1 * x - y
    cost = r' * r
    g = AS1' * r
    cost_change = Inf

    # ----- Initialization loop ----- #

    while !isempty(free_set)

        AS1_free = AS1[:, free_set]
        y_free = y - AS1 * (x .* active_set)
        z = AS1_free \ y_free
        lbv = (z .< lb[free_set])
        ubv = (z .> ub[free_set])
        v = lbv .| ubv

        if any(lbv)
            ind = free_set[lbv]
            x[ind] = lb[ind]
            active_set[ind] .= true
            on_bound[ind] .= -1
        end
        if any(ubv)
            ind = free_set[ubv]
            x[ind] = ub[ind]
            active_set[ind] .= true
            on_bound[ind] .= 1
        end

        ind = free_set[.!v]
        x[ind] = z[.!v]
        r = AS1 * x - y
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
            x_free = x[free_set]
            lb_free = lb[free_set]
            ub_free = ub[free_set]
            AS1_free = AS1[:,free_set]
            y_free = y - AS1 * (x .* active_set)
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
                x[free_set] = x_free

                if i < length(lbv)
                    on_bound[free_set[i_free]] = -1
                else
                    on_bound[free_set[i_free]] = 1
                end
            else
                x_free = z
                x[free_set] = x_free
                break
            end
        end

        r = AS1 * x - y
        cost_new = 0.5 * (r' * r)
        cost_change = cost - cost_new

        (cost_change < bnbparams.tol * cost) && break

        cost = cost_new
        g = AS1' * r
        optimality = compute_kkt_optimality(g,on_bound)
    end

    obj = 0.5 * norm(r,2)^2 + λ * n
    
    return obj, x
end
