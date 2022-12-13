module BnbScreening

using Dates
using Printf
using Random
using CPLEX
using Distributions
using JuMP
using LinearAlgebra

version() = "v1.0"
author() = "Theo Guyard"
contact() = "theo.guyard@insa-rennes.fr"
license() = "AGPL 3.0 License (https://www.gnu.org/licenses/agpl-3.0.en.html)"

# Module files
include("types.jl")
include("activeset.jl")
include("bnb.jl")
include("bvls.jl")
include("data.jl")
include("displays.jl")
include("mips.jl")
include("screen.jl")

export data_gaussian, data_toeplitz
export LbMethod, UbMethod
export solve_bnb, BnbParams, BnbResults
export solve_mip, MipParams, MipResults

end
