mutable struct DDOptimizer <: AbstractGraphOptimizer
    dd_model::DDModel
end

JuMP.optimize!(graph::ModelGraph,optimizer::DDOptimizer) = dual_decomposition_solve(graph)
