using Revise
using JuMP
using GLPK
using ModelGraphs

include("DualDecompositionSolver.jl")
using .DualDecompositionSolver

const DDSolver = DualDecompositionSolver

m2 = Model(with_optimizer(GLPK.Optimizer))
@variable(m2, xs[i in 1:2], Bin)
@variable(m2, y[i in 1:2], Bin)
@constraint(m2, y[1] + y[2] <= 1)
@constraint(m2, 8xs[1] + 2xs[2] + y[1] + 4y[2] <= 10)
@objective(m2, Max, 4y[2])

m3 = Model(with_optimizer(GLPK.Optimizer))
@variable(m3, xs[i in 1:2], Bin)
@variable(m3, y[i in 1:2], Bin)
@constraint(m3, y[1] + y[2] <= 1)
@constraint(m3, 8xs[1] + 2xs[2] + y[1] + 4y[2] <= 10)
@objective(m3, Max, y[2])

## Model Graph
graph = ModelGraph()
JuMP.set_optimizer(graph.mastermodel,with_optimizer(GLPK.Optimizer))

@linkvariable(graph, xm[i in 1:2],Bin)
@masterconstraint(graph, xm[1] + xm[2] <= 1)
@objective(graph.mastermodel, Max, 16xm[1].vref + 10xm[2].vref)

n1 = add_node!(graph,m2)
n2 = add_node!(graph,m3)

for i = 1:2
    link_variables!(graph[:xm][i],n1[:xs][i])
    link_variables!(graph[:xm][i],n2[:xs][i])
end

# link constraints between models
# @linkconstraint(graph, [i in 1:2], n1[:xm][i] == n2[:xs][i])
# @linkconstraint(graph, n3[:x3][1] + n1[:xm][1] + n2[:xs][1] <= 2)
# @linkconstraint(graph, n1[:xm][2] <= n3[:y][1])

# dd_model = DDModel(graph)
# DDSolver.optimize!(dd_model)
#
# #print solution
# println("Dual decomposition solution: ")
# solution = Float64[]
# for subproblem in dd_model.subproblems
#     append!(solution,value.(all_variables(subproblem)))
# end
# println(solution)
#
# #verify with full problem
glpk = with_optimizer(GLPK.Optimizer)
m,ref_map = aggregate(graph)
optimize!(m,glpk)

println()
println("Pure glpk solution:")
println(JuMP.value.(all_variables(m)))
