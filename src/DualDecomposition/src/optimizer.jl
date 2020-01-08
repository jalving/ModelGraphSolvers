mutable struct DDOptimizer <: AbstractGraphOptimizer
    args::Vector
    kwargs::Dict
    ddmodel::Union{Nothing,DDModel}
end

function DDOptimizer(args...;kwargs...)
    return DDOptimizer(args,kwargs,nothing)
end

JuMP.optimize!(graph::ModelGraph,optimizer::DDOptimizer) = dual_decomposition_solve(graph,optimizer.args...;optimizer.kwargs...)
