using Distributed
if nprocs() == 1
    addprocs(2)
end
@everywhere using Pkg
@everywhere using Revise
@everywhere Pkg.activate(".")
@everywhere using ModelGraphs
@everywhere using JuMP
@everywhere using GLPK

include("model_link_constraint.jl")

dd_model = DDModel(graph)
# for sub in dd_model.subproblems
#     JuMP.set_optimizer(sub,with_optimizer(GLPK.Optimizer))
# end

parallel_solve(dd_model,workers())

#print solution
println("Distirbuted Dual decomposition solution: ")


parallel_solve(dd_model)


# # #verify with full problem
# glpk = with_optimizer(GLPK.Optimizer)
# m,ref_map = aggregate(graph)
# optimize!(m,glpk)
#
# println()
# println("Pure glpk solution:")
# println(JuMP.value.(all_variables(m)))


solution = Float64[]
for subproblem in dd_model.subproblems
    append!(solution,value.(all_variables(subproblem)))
end
println(solution)
