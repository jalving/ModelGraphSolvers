module DDSolver

using AlgebraicGraphs
import AlgebraicGraphs.AbstractGraphOptimizer
using JuMP

using MathOptInterface
const MOI = MathOptInterface

using SparseArrays
using LinearAlgebra

export DDModel, dual_decomposition_solve

include("utils.jl")

include("solution.jl")

include("lagrange_model.jl")

include("multiplier_updates.jl")

end
