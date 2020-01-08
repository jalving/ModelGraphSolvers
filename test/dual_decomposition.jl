using JuMP
using GLPK
using ModelGraphs
using ModelGraphSolvers

m1 = Model(with_optimizer(GLPK.Optimizer))

@variable(m1, xm[i in 1:2],Bin)
@constraint(m1, xm[1] + xm[2] <= 1)
@objective(m1, Max, 16xm[1] + 10xm[2])

m2 = Model(with_optimizer(GLPK.Optimizer))
@variable(m2, xs[i in 1:2],Bin)
@variable(m2, y[i in 1:2], Bin)
@constraint(m2, y[1] + y[2] <= 1)
@constraint(m2, 8xs[1] + 2xs[2] + y[1] + 4y[2] <= 10)
@objective(m2, Max, 4y[2])

## Model Graph
graph = ModelGraph()

n1 = add_node!(graph)
set_model(n1,m1)
n2 = add_node!(graph)
set_model(n2,m2)

## Linking
# m1[x] = m2[x]  ∀i ∈ {1,2}
@linkconstraint(graph, [i in 1:2], n1[:xm][i] == n2[:xs][i])

dual_decomposition_solve(graph)

true
