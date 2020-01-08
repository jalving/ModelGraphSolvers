module DualDecomposition

using ModelGraphs
using JuMP
using MathOptInterface
const MOI = MathOptInterface
using GLPK

using SparseArrays
using LinearAlgebra
using Distributed

export DDModel, dual_decomposition_solve, DDOptimizer

include("../../utils.jl")

#include("../../distribute.jl")

include("solution.jl")

include("dual_decomp_model.jl")

include("multiplier_updates.jl")

include("optimizer.jl")

end
