using Revise
using JuMP
using GLPK
using ModelGraphs

include("DualDecompositionSolver.jl")
using .DualDecompositionSolver

DDSolver = DualDecompositionSolver

m2 = Model(with_optimizer(GLPK.Optimizer))
@variable(m2, xs[i in 1:2])
@variable(m2, y[i in 1:2], Bin)
@constraint(m2, y[1] + y[2] <= 1)
@constraint(m2, xs[1] + xs[2] + y[1] + 4y[2] <= 10)
@objective(m2, Min, xs[1] + xs[2])

m3 = Model(with_optimizer(GLPK.Optimizer))
@variable(m3, xs[i in 1:2])
@variable(m3, y[i in 1:2], Bin)
@constraint(m3, y[1] + y[2] <= 1)
@constraint(m3, xs[1] + xs[2] + y[1] + 4y[2] <= 10)
@objective(m3, Min, y[2])

## Model Graph
graph = ModelGraph()
JuMP.set_optimizer(graph.mastermodel,with_optimizer(GLPK.Optimizer))

@linkvariable(graph, xm[i in 1:2] >= 2, Int)
@masterconstraint(graph, xm[1] + xm[2] <= 10)
#@objective(graph.mastermodel, Max, 16*xm[1].vref + 10*xm[2].vref)

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

#
# #verify with full problem
glpk = with_optimizer(GLPK.Optimizer)
m,ref_map = aggregate(graph)
optimize!(m,glpk)

println()
println("Pure glpk solution:")
println(JuMP.value.(all_variables(m)))


dd_model = DDSolver.DDModel(graph)

for subproblem in dd_model.subproblems
    JuMP.set_optimizer(subproblem,glpk)
end

DDSolver.optimize!(dd_model)
#print solution
println("Dual decomposition solution: ")
for subproblem in dd_model.subproblems
    println(value.(all_variables(subproblem)))
end


# link_eq_matrix,b_eq,link_eq_variables = DDSolver.prepare_link_matrix(link_eq_constraints,sub_var_map)
# link_ineq_matrix,b_ineq,link_ineq_variables, = prepare_link_matrix(link_ineq_constraints,sub_var_map)
# link_var_matrix,link_var_variables = prepare_link_matrix(link_variables,master_var_map)

#link_var_matrix,link_var_variables,link_var_map = DDSolver.prepare_link_matrix(linkvariables)
