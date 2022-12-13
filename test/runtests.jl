using BnbScreening
using Test

using MathOptInterface

const MOI = MathOptInterface

@testset "gaussian" begin
    k, m, n, snr = 5, 500, 1000, 10*log10(10)
    A, y, λ, M = data_gaussian(k,m,n,snr)
    result_mip = solve_mip(A, y, λ, M)
    result_bnb = solve_bnb(A, y, λ, M, bnbparams=BnbParams(verbosity=false))
    result_bnb_scr = solve_bnb(A, y, λ, M, bnbparams=BnbParams(verbosity=false,screening=true))
    @test result_mip.termination_status == MOI.OPTIMAL
    @test result_bnb.termination_status == MOI.OPTIMAL
    @test result_bnb_scr.termination_status == MOI.OPTIMAL
    @test abs(result_mip.objective_value - result_bnb.objective_value) < 1e-6
    @test abs(result_mip.objective_value - result_bnb_scr.objective_value) < 1e-6
    @test abs(result_bnb.objective_value - result_bnb_scr.objective_value) < 1e-6
end

@testset "toeplitz" begin
    k, m, n, snr = 5, 500, 300, 10*log10(10)
    A, y, λ, M = data_toeplitz(k,m,n,snr)
    result_mip = solve_mip(A, y, λ, M)
    result_bnb = solve_bnb(A, y, λ, M, bnbparams=BnbParams(verbosity=false))
    result_bnb_scr = solve_bnb(A, y, λ, M, bnbparams=BnbParams(verbosity=false,screening=true))
    @test result_mip.termination_status == MOI.OPTIMAL
    @test result_bnb.termination_status == MOI.OPTIMAL
    @test result_bnb_scr.termination_status == MOI.OPTIMAL
    @test abs(result_mip.objective_value - result_bnb.objective_value) < 1e-6
    @test abs(result_mip.objective_value - result_bnb_scr.objective_value) < 1e-6
    @test abs(result_bnb.objective_value - result_bnb_scr.objective_value) < 1e-6
end
