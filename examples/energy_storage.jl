using JuMP, Ipopt
using ModelGraphs

# Function to create individual node models
function storage(M,d)
    m = Model(with_optimizer(Ipopt.Optimizer),print_level=0))
    @variable(m, x[1:M+1])
    @variable(m, u[1:M])
    @constraint(m, dynamics[i in 1:M], x[i+1] == x[i] + u[i] + d[i])
    @objective(m, Min, 0.001*sum(x[i]^2 for i in 1:M) + sum(u[i]^2 for i in 1:M))
    return m
end

n_nodes = 50
N = 800
k = round(Int64,n_nodes/N)
d = sin.(1:N)

disturbance_matrix = reshape(d, (k, div(length(d), k)))
graph = ModelGraph()

for i = 1:n_nodes
    node_disturbance = disturbance_matrix[:,i]
    add_node!(graph,storagemodel(N,node_disturbance))
end

n1 = getnode(graph,1)
@constraint(n1,n1[:x][1] == 0)                      #First node in planning horizon has initial condition
for i = 1:n_nodes - 1
    ni = getnode(graph,i)
    nj = getnode(graph,i+1)
    @linkconstraint(graph, ni[:x][M] == nj[:x][1])  #last state in partition i is first state in partition j
end
