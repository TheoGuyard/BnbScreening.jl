tolsign(x::Float64, tol::Float64) = sign(x) * (abs(x) > tol)

function symbi_app(Ainv::Matrix{Float64}, a::Vector{Float64})
    Ainv_updated = Matrix{Float64}(undef,length(a),length(a))
    @views begin
        d = a[1:end-1]
        g = a[end]
        dw = BLAS.symv('U',Ainv,-d)
        gw = 1. / (g + dw' * d)
        gw_dw = gw .* dw
        BLAS.syr!('U',gw,dw,Ainv)
        Ainv_updated[1:end-1,1:end-1] = Ainv
        Ainv_updated[1:end-1,end] = gw_dw
        Ainv_updated[end,end] = gw
    end
    return Ainv_updated
end

function symbi_rem(Ainv::Matrix{Float64}, i::Int)
    n = size(Ainv)[1]
    Ainv_updated = Matrix{Float64}(undef,n-1,n-1)
    @views begin
        q = -1. / Ainv[i,i]
        c1 = 1:i-1
        r1 = Ainv[c1,i]
        c2 = i+1:n
        r2 = Ainv[i,c2]
        c2red = i:n-1
        Ainv_updated[c1,c1] = BLAS.syr!('U',q,r1,Ainv[c1,c1])
        Ainv_updated[c1,c2red] = BLAS.ger!(q,r1,r2,Ainv[c1,c2])
        Ainv_updated[c2red,c2red] = BLAS.syr!('U',q,r2,Ainv[c2,c2])
    end
    return Ainv_updated
end

function symbi_ro!(Ainv::Matrix{Float64}, r::Vector{Float64}, g::Float64)
    Ainv_r = BLAS.symv('U',Ainv,r)
    c = -g / (1. + g * (r' * Ainv_r))
    BLAS.syr!('U',c,Ainv_r,Ainv)
    return nothing
end

function swap_Sbarin_to_Sbar0!(
    node::Node, 
    asviews::ActivesetViews, 
    A::Matrix{Float64}, 
    i::Int
    )
    iSbarin = findfirst(j->j==i,node.support.Sbarin)
    node.GSbarin_inv = symbi_rem(node.GSbarin_inv,iSbarin)
    push!(node.support.Sbar0,i)
    deleteat!(node.support.Sbarin,iSbarin)
    deleteat!(asviews.xSbarin,iSbarin)
    deleteat!(asviews.sSbarin,iSbarin)
    deleteat!(asviews.xnewSbarin,iSbarin)
    asviews.ASbarin = A[:,node.support.Sbarin]
    return nothing
end

function swap_Sbarin_to_Sbarbox!(node::Node, 
    asviews::ActivesetViews, 
    A::Matrix{Float64}, 
    i::Int
    )
    iSbarin = findfirst(j->j==i,node.support.Sbarin)
    node.GSbarin_inv = symbi_rem(node.GSbarin_inv,iSbarin)
    push!(node.support.Sbarbox,i)
    push!(asviews.xSbarbox,popat!(asviews.xSbarin,iSbarin))
    push!(asviews.sSbarbox,popat!(asviews.sSbarin,iSbarin))
    deleteat!(asviews.xnewSbarin,iSbarin)
    deleteat!(node.support.Sbarin,iSbarin)
    asviews.ASbarin = A[:,node.support.Sbarin]
    asviews.ASbarbox = A[:,node.support.Sbarbox]
    return nothing
end

function swap_S1in_to_S1box!(
    node::Node, 
    asviews::ActivesetViews, 
    A::Matrix{Float64}, 
    B::Matrix{Float64}, 
    i::Int
    )
    iS1in = findfirst(j->j==i,node.support.S1in)
    push!(node.support.S1box,i)
    push!(asviews.xS1box,popat!(asviews.xS1in,iS1in))
    push!(asviews.sS1box,popat!(asviews.sS1in,iS1in))
    deleteat!(asviews.xnewS1in,iS1in)
    deleteat!(node.support.S1in,iS1in)
    asviews.AS1in = A[:,node.support.S1in]
    asviews.AS1box = A[:,node.support.S1box]

    @views begin
        a = A[:,i]
        d = vcat(node.GS1in_inv[1:(iS1in-1),iS1in], node.GS1in_inv[iS1in,(iS1in+1):end])
        e = node.GS1in_inv[iS1in,iS1in]
        Ad = asviews.AS1in * d
        r1 = (Ad./e) + a
        BLAS.syr!('U',e,r1,B)
        node.GS1in_inv = symbi_rem(node.GS1in_inv,iS1in)
        symbi_ro!(node.GSbarin_inv, asviews.ASbarin'*r1, e)
    end
    
    return nothing
end

function swap_Sbar0_to_Sbarin!(
    node::Node, 
    asviews::ActivesetViews, 
    A::Matrix{Float64}, 
    B::Matrix{Float64}, 
    i::Int
    )
    @views begin  
        a = A[:,i]
        Ba = BLAS.symv('U',B,a)
        ABa = asviews.ASbarin' * Ba
        aBa = a' * Ba
        v = vcat(ABa,aBa)
        node.GSbarin_inv = symbi_app(node.GSbarin_inv,v)
    end
    push!(node.support.Sbarin,i)
    filter!(j->j!=i,node.support.Sbar0)
    push!(asviews.xSbarin, 0.)
    push!(asviews.sSbarin, 0.)
    push!(asviews.xnewSbarin, 0.)
    asviews.ASbarin = A[:,node.support.Sbarin]
    return nothing
end

function swap_Sbarbox_to_Sbarin!(
    node::Node, 
    asviews::ActivesetViews, 
    A::Matrix{Float64}, 
    B::Matrix{Float64}, 
    i::Int
    )
    iSbarbox = findfirst(j->j==i,node.support.Sbarbox)
    @views begin
        a = A[:,i]
        Ba = BLAS.symv('U',B,a)
        ABa = asviews.ASbarin' * Ba
        aBa = a' * Ba
        v = vcat(ABa,aBa)
        node.GSbarin_inv = symbi_app(node.GSbarin_inv,v)
    end
    push!(node.support.Sbarin, popat!(node.support.Sbarbox,iSbarbox))
    push!(asviews.xSbarin, popat!(asviews.xSbarbox,iSbarbox))
    push!(asviews.sSbarin, popat!(asviews.sSbarbox,iSbarbox))
    push!(asviews.xnewSbarin, asviews.xSbarin[end])
    asviews.ASbarin = A[:,node.support.Sbarin]
    asviews.ASbarbox = A[:,node.support.Sbarbox]
    return nothing
end

function swap_S1box_to_Sbarin!(
    node::Node, 
    asviews::ActivesetViews, 
    A::Matrix{Float64}, 
    B::Matrix{Float64}, 
    i::Int
    )
    iS1box = findfirst(j->j==i,node.support.S1box)
    @views begin
        a = A[:,i]
        Ba = BLAS.symv('U',B,a)
        ABa = asviews.ASbarin' * Ba
        aBa = a' * Ba
        v = vcat(ABa,aBa)
        node.GSbarin_inv = symbi_app(node.GSbarin_inv,v)
    end
    push!(node.support.Sbarin, popat!(node.support.S1box,iS1box))
    push!(asviews.xSbarin, popat!(asviews.xS1box,iS1box))
    push!(asviews.sSbarin, popat!(asviews.sS1box,iS1box))
    push!(asviews.xnewSbarin, asviews.xSbarin[end])
    asviews.ASbarin = A[:,node.support.Sbarin]
    asviews.AS1box = A[:,node.support.S1box]
    return nothing
end

function swap_Sbar_to_S0!(
    node::Node, 
    asviews::ActivesetViews, 
    A::Matrix{Float64}, 
    i::Int
    )
    iSbarin = findfirst(j->j==i,node.support.Sbarin)
    if !isa(iSbarin,Nothing)
        node.GSbarin_inv = symbi_rem(node.GSbarin_inv,iSbarin)
        deleteat!(node.support.Sbarin,iSbarin)
        deleteat!(asviews.xSbarin,iSbarin)
        deleteat!(asviews.sSbarin,iSbarin)
        deleteat!(asviews.xnewSbarin,iSbarin)
        asviews.ASbarin = A[:,node.support.Sbarin]
    else
        iSbarbox = findfirst(j->j==i,node.support.Sbarbox)
        if !isa(iSbarbox,Nothing)
            deleteat!(node.support.Sbarbox,iSbarbox)
            deleteat!(asviews.xSbarbox,iSbarbox)
            deleteat!(asviews.sSbarbox,iSbarbox)
            asviews.ASbarbox = A[:,node.support.Sbarbox]
        else
            filter!(j->j!=i,node.support.Sbar0)
        end
    end
    return nothing
end

function swap_Sbar_to_S1!(
    node::Node, 
    asviews::ActivesetViews, 
    A::Matrix{Float64}, 
    B::Matrix{Float64}, 
    i::Int
    )
    iSbarbox = findfirst(j->j==i,node.support.Sbarbox)
    if !isa(iSbarbox,Nothing)
        push!(node.support.S1box, popat!(node.support.Sbarbox,iSbarbox))
        push!(asviews.xS1box, popat!(asviews.xSbarbox,iSbarbox))
        push!(asviews.sS1box, popat!(asviews.sSbarbox,iSbarbox))
        asviews.ASbarbox = A[:,node.support.Sbarbox]
        asviews.AS1box = A[:,node.support.S1box]
    else
        @views begin
            asviews.AS1in = A[:,node.support.S1in]
            a = A[:,i]
            Ata = asviews.AS1in' * a
            ata = a' * a
            d = - BLAS.symv('U',node.GS1in_inv,Ata)
            g = -1. / (ata + d'*Ata)
            r1 = asviews.AS1in * d + a
            v = vcat(Ata,ata)
            node.GS1in_inv = symbi_app(node.GS1in_inv,v)
            BLAS.syr!('U',g,r1,B)
        end
        iSbarin = findfirst(j->j==i,node.support.Sbarin)
        if !isa(iSbarin,Nothing)
            node.GSbarin_inv = symbi_rem(node.GSbarin_inv,iSbarin)
            push!(node.support.S1in, popat!(node.support.Sbarin,iSbarin))
            push!(asviews.xS1in, popat!(asviews.xSbarin,iSbarin))
            push!(asviews.sS1in, popat!(asviews.sSbarin,iSbarin))
            push!(asviews.xnewS1in, popat!(asviews.xnewSbarin,iSbarin))
            asviews.ASbarin = A[:,node.support.Sbarin]
            asviews.AS1in = A[:,node.support.S1in]
            symbi_ro!(node.GSbarin_inv, asviews.ASbarin'*r1,g)
        else
            iSbar0 = findfirst(j->j==i,node.support.Sbar0)
            push!(node.support.S1in, popat!(node.support.Sbar0,iSbar0))
            push!(asviews.xS1in, 0.)
            push!(asviews.sS1in, 0.)
            push!(asviews.xnewS1in, 0.)
            asviews.AS1in = A[:,node.support.S1in]
            symbi_ro!(node.GSbarin_inv, asviews.ASbarin'*r1,g)
        end
    end
    return nothing
end

function solve_kkt_system!(
    node::Node,
    asviews::ActivesetViews,
    B::Matrix{Float64},
    y::Vector{Float64},
    λ::Float64,
    M::Float64
    )
    Axbox = asviews.ASbarbox * asviews.xSbarbox + asviews.AS1box * asviews.xS1box
    r = y - Axbox
    rhsSbarin = asviews.ASbarin' * BLAS.symv('U',B,r) - (λ/M) .* asviews.sSbarin
    asviews.xnewSbarin = BLAS.symv('U',node.GSbarin_inv,rhsSbarin)
    rhsS1in = asviews.AS1in' * (r - asviews.ASbarin * asviews.xnewSbarin)
    asviews.xnewS1in = BLAS.symv('U',node.GS1in_inv,rhsS1in)
    return Axbox
end

function discrete_line_search!(
    node::Node,
    asviews::ActivesetViews,
    A::Matrix{Float64},
    Axbox::Vector{Float64},
    y::Vector{Float64},
    λ::Float64,
    M::Float64,
    bnbparams::BnbParams,
    )

    # TODO : optimize
    
    xdiffSbarin = asviews.xnewSbarin - asviews.xSbarin
    xdiffS1in = asviews.xnewS1in - asviews.xS1in
    Ax = asviews.ASbarin * asviews.xSbarin + asviews.AS1in * asviews.xS1in + Axbox
    Axdiff = asviews.ASbarin * xdiffSbarin + asviews.AS1in * xdiffS1in
    u = y - Ax
    λ_M = λ / M
    λ_S1 = λ * length(node.branch.S1)

    a = Axdiff' * Axdiff
    b_base = u' * Axdiff

    if a == 0.
        Atu = A' * u
        lb = 0.5*u'*u + λ_M*asviews.xSbarin'*asviews.sSbarin + λ_S1
        return u, Atu, lb
    end

    tflip = []
    tmax  = 1.
    for i in 1:length(node.support.Sbarin)
        if xdiffSbarin[i] < 0.
            tmax = min(tmax, -(M+asviews.xSbarin[i])/xdiffSbarin[i])
        else
            tmax = min(tmax, (M-asviews.xSbarin[i])/xdiffSbarin[i])
        end
        push!(tflip, -asviews.xSbarin[i]/xdiffSbarin[i])
    end
    for i in 1:length(node.support.S1in)
        if xdiffS1in[i] < 0.
            tmax = min(tmax, -(M+asviews.xS1in[i])/xdiffS1in[i])
        else
            tmax = min(tmax, (M-asviews.xS1in[i])/xdiffS1in[i])
        end
        push!(tflip, -asviews.xS1in[i]/xdiffS1in[i])
    end

    filter!(ti -> 0. < ti < tmax, sort!(unique!(tflip)))
    pushfirst!(tflip,0.)
    push!(tflip,tmax)

    topt = 0.
    uopt = u
    vopt = Inf
    
    for i in 1:(length(tflip)-1)

        tl   = tflip[i]
        tu   = tflip[i+1]
        tmed = (tl+tu)/2.

        thetaSbarin = sign.(asviews.xSbarin + tmed .* xdiffSbarin)

        b = b_base - λ_M * thetaSbarin' * xdiffSbarin

        tinf = max(tl,min(tu,b/a))
        uinf = u - tinf .* Axdiff
        vinf = 0.5 * norm(uinf,2)^2 + λ_M * norm(asviews.xSbarin + tinf .* xdiffSbarin,1)

        if vinf < vopt
            topt = tinf
            uopt = uinf
            vopt = vinf
        end
    end

    asviews.xSbarin += topt .* xdiffSbarin
    asviews.xS1in += topt .* xdiffS1in
    asviews.sSbarin = tolsign.(asviews.xSbarin,bnbparams.tol)
    asviews.sS1in = tolsign.(asviews.xS1in,bnbparams.tol)
    Atuopt = A' * uopt
    # lb = vopt + (λ_S1 + λ * sum(node.support.Sbarbox))
    w = y - u
    pSb = @. (M * abs(Atuopt[node.branch.Sbar]) - λ)
    pS1 = @. (M * abs(Atuopt[node.branch.S1]) - λ)
    lb = 0.5 * (y' * y - w' * w) - sum(max.(pSb, 0.)) - sum(pS1)

    return uopt, Atuopt, lb
end

function detect_blocking_behaviour(swapped::Vector{Int}, imax::Int)
    return (length(swapped)==1) & (imax in swapped)
end

function active_optimality(
    asviews::ActivesetViews,
    M::Float64,
    bnbparams::BnbParams
    )
    cond = (
        all(asviews.sSbarin .== tolsign.(asviews.xnewSbarin,bnbparams.tol)) & 
        all(asviews.sS1in .== tolsign.(asviews.xnewS1in,bnbparams.tol)) &
        (norm(asviews.xnewSbarin,Inf) - M <= bnbparams.tol) &
        (norm(asviews.xnewS1in,Inf) - M <= bnbparams.tol) 
    )
    return cond 
end

function inactive_optimality(
    node::Node,
    asviews::ActivesetViews,
    Atu::Vector{Float64},
    λ::Float64,
    M::Float64,
    bnbparams::BnbParams
    )
    vmax = 0.
    imax = nothing
    smax = ""
    λ_M = λ / M
    λ_M_m_tol =  λ_M - bnbparams.tol
    λ_M_p_tol =  λ_M + bnbparams.tol
    m_tol = -bnbparams.tol
    for i in node.support.Sbar0
        v = abs(Atu[i]) - λ_M_p_tol
        if v > vmax
            vmax = v
            imax = i
            smax = "Sbar0"
        end
    end
    for (iSbarbox,i) in enumerate(node.support.Sbarbox)
        v = λ_M_m_tol - Atu[i] * asviews.sSbarbox[iSbarbox]
        if v > vmax
            vmax = v
            imax = i
            smax = "Sbarbox"
        end
    end
    for (iS1box,i) in enumerate(node.support.S1box)
        v = m_tol - Atu[i] * asviews.sS1box[iS1box]
        if v > vmax
            vmax = v
            imax = i
            smax = "S1box"
        end
    end
    return imax, smax
end

function new_tentative_sign!(
    node::Node, 
    asviews::ActivesetViews, 
    imax::Int, 
    smax::String, 
    A::Matrix{Float64}, 
    B::Matrix{Float64}
    )
    if smax == "Sbar0"
        swap_Sbar0_to_Sbarin!(node,asviews,A,B,imax)
    elseif smax == "Sbarbox"
        swap_Sbarbox_to_Sbarin!(node,asviews,A,B,imax)
    elseif smax == "S1box"
        swap_S1box_to_Sbarin!(node,asviews,A,B,imax)
    end
    return nothing
end

function update_sets!(
    node::Node,
    asviews::ActivesetViews, 
    A::Matrix{Float64}, 
    B::Matrix{Float64},
    M::Float64,
    bnbparams::BnbParams
    )

    swapped = Vector{Int}()
    for (iS1in,i) in enumerate(node.support.S1in)
        absxi = abs(asviews.xS1in[iS1in])
        if abs(M - absxi) <= bnbparams.tol
            swap_S1in_to_S1box!(node,asviews,A,B,i)
            push!(swapped,i)
        end
    end
    for (iSbarin,i) in enumerate(node.support.Sbarin)
        absxi = abs(asviews.xSbarin[iSbarin])
        if absxi <= bnbparams.tol
            swap_Sbarin_to_Sbar0!(node,asviews,A,i)
            push!(swapped,i)
        elseif abs(M - absxi) <= bnbparams.tol
            swap_Sbarin_to_Sbarbox!(node,asviews,A,i)
            push!(swapped,i)
        end
    end
    return swapped
end

function solve_relax_activeset!(
    tree::Tree,
    node::Node,
    A::Matrix{Float64},
    y::Vector{Float64},
    λ::Float64,
    M::Float64,
    bnbparams::BnbParams,
    )

    asviews = ActivesetViews(node,A,bnbparams)
    B = I(size(A)[1]) - BLAS.symm('R','U',node.GS1in_inv,asviews.AS1in) * asviews.AS1in'
    node.branch_on[2] == 0 && swap_Sbar_to_S0!(node,asviews,A,node.branch_on[1])
    node.branch_on[2] == 1 && swap_Sbar_to_S1!(node,asviews,A,B,node.branch_on[1])

    it = 0
    lb = Inf

    while true
        it += 1
        it >= bnbparams.maxiter && error("activetset : iteration limit")

        Axbox = solve_kkt_system!(node,asviews,B,y,λ,M)
        u, Atu, lb = discrete_line_search!(node,asviews,A,Axbox,y,λ,M,bnbparams)
        swapped = update_sets!(node,asviews,A,B,M,bnbparams)

        if bnbparams.screening
            prune = screen!(tree,node,asviews,u,Atu,B,A,y,λ,M)
            prune && break
        end

        if active_optimality(asviews,M,bnbparams)
            imax, smax = inactive_optimality(node,asviews,Atu,λ,M,bnbparams)
            isa(imax,Nothing) && break
            new_tentative_sign!(node,asviews,imax,smax,A,B)
            detect_blocking_behaviour(swapped,imax) && break
        end
    end

    node.lb = lb
    node.x[node.support.Sbar0] .= 0.
    node.x[node.support.Sbarin] = asviews.xSbarin
    node.x[node.support.Sbarbox] = asviews.xSbarbox
    node.x[node.support.S1in] = asviews.xS1in
    node.x[node.support.S1box] = asviews.xS1box
    
    return nothing
end
