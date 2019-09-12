using ModelGraphs
using JuMP
using GLPK

graph = ModelGraph()

function sub_model()
    m1 = Model()
    @variable(m1,x[1:2] >= 0)
    @constraint(m1,x[1] >= 3)
    @objective(m1,Min,x[1] + x[2])
    setsolver(m1,GLPKSolverLP())
    return m1
end

n2 = add_node!(graph,sub_model())
n3 = add_node!(graph,sub_model())

@linkvariable(graph,x[1:2] >= 0)
link_variables!(graph[:x][1],n2[:x][1])
link_variables!(graph[:x][2],n2[:x][2])
link_variables!(graph[:x][2],n3[:x][2])


optimize!(graph,with_optimizer(BendersOptimizer))
