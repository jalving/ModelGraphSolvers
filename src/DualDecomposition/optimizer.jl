mutable struct DDOptimizer <: AbstractGraphOptimizer
    dd_model::DDModel
end

function DDOptimizer()

end

JuMP.optimize!(graph::ModelGraph,optimizer::DDOptimizer) = dual_decomposition_solve(graph,optimizer.args...;optimizer.kwargs...)
