# using .DualDecomposition
using Pkg
Pkg.activate(".")
include("../src/DualDecomposition.jl")

using Revise
using JuMP
using GLPK
using ModelGraphs

include("model_link_constraint.jl")

dd_model = DDModel(graph)
for sub in dd_model.subproblems
    JuMP.set_optimizer(sub,with_optimizer(GLPK.Optimizer))
end

solve(dd_model)

#print solution
println("Dual decomposition solution: ")
solution = Float64[]
for subproblem in dd_model.subproblems
    append!(solution,value.(all_variables(subproblem)))
end
println(solution)

# #verify with full problem
glpk = with_optimizer(GLPK.Optimizer)
m,ref_map = aggregate(graph)
optimize!(m,glpk)

println()
println("Pure glpk solution:")
println(JuMP.value.(all_variables(m)))


#dual_decomposition_solve(graph)
