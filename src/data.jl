function validate_data(
    A::Array{Float64,2}, 
    y::Array{Float64,1}, 
    λ::Float64, 
    M::Float64
    )
    size(A)[1] == size(y)[1] || throw(ArgumentError("Dimension missmatch"))
    λ > 0 || throw(ArgumentError("λ must be positive"))
    M > 0 || throw(ArgumentError("M must be positive"))
end

# ============================================================================ #

"""
    data_gaussian(
        k::Int=5, 
        m::Int=500, 
        n::Int=1000, 
        snr::Float64=10.
        )

Generate a Gaussian setup.
"""
function data_gaussian(
    k::Int=5, 
    m::Int=500, 
    n::Int=1000, 
    snr::Float64=10.
    )

    x = zeros(n)
    a = rand(Normal(0,1),k)
    a += sign.(a)
    x[randperm(n)[1:k]] = a

    A = rand(Normal(0,1),m,n)
    for i in 1:n
        A[:,i] /= norm(A[:,i],2)
    end

    σ = norm(A*x,2) / sqrt(m * 10^(snr/10))
    ϵ = rand(Normal(0,σ),m)
    y = A * x .+ ϵ
    λ = 2 * σ * log(n/k - 1)
    M = 1.5 * norm(A'*y,Inf)

    return A, y, λ, M
end

"""
    data_gaussian(
        k::Int=5, 
        m::Int=500, 
        n::Int=1000, 
        snr::Float64=10.
        )

Generate a Toeplitz setup.
"""
function data_toeplitz(
    k::Int=5, 
    m::Int=500, 
    n::Int=300, 
    snr::Float64=10.
    )

    x = zeros(n)
    a = rand(Normal(0,1),k)
    a += sign.(a)
    x[randperm(n)[1:k]] = a

    tmax = n/2
    t = collect(1:m)
    A = zeros(m,n)
    offset = (n-m)/2
    for i in 1:n
        A[:,i] = @. sinc((offset + (t - i)) * (tmax/m))
        A[:,i] /= norm(A[:,i],2)
    end
    A += rand(m,n)./10000

    σ = norm(A*x,2) / sqrt(m * 10^(snr/10))
    ϵ = rand(Normal(0,σ),m)
    y = A * x .+ ϵ
    λ = 2 * σ * log(n/k - 1)
    M = 1.5 * norm(A'*y,Inf)

    return A, y, λ, M
end
