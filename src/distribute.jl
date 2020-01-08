#Distribute subproblems on a DD Model to workers.  return references to subproblems
function distribute_subproblems(subproblems::Vector{JuMP.Model},to_workers::Vector{Int64})

    println("Distributing $(length.dd_model.subproblems) among $(length(to_workers))...")

    n_problems = length(subproblems)
    n_workers = length(to_workers)
    problems_per_worker = Int64(floor(n_problems/n_workers))

    #Allocate modelnodes onto provided workers
    allocations = []
    j = 1
    while  j <= n_problems
        if j + problems_per_worker >= n_problems
            push!(allocations,subproblems[j:end])
            break
        else
            push!(allocations,subproblems[j:j+problems_per_worker - 1])
        end
        global j += problems_per_worker
    end
    worker_map = Dict(zip(to_workers,allocations))

    #Fill channels with subproblems
    subproblem_map = Dict()
    worker_channel = Dict()
    for (worker,allocation) in worker_map
        channel = RemoteChannel(1)
        worker_channel[worker] = channel
        for subproblem in allocation
            subproblem_map[subproblem] = worker
            @spawnat(1, put!(channel, subproblem))
        end
    end

    remote_subproblem_references = []
    remote_mapping = Dict()

    @sync begin
        for (subproblem,worker) in subproblem_map
            channel = worker_channel[worker]
            sub_reference = @spawnat worker Core.eval(Main, Expr(:(=), gensym(), take!(channel)))   #wait(ref1)
            push!(remote_subproblem_references,sub_reference)
            remote_mapping[sub_reference] = worker
        end
    end

    return remote_subproblem_references, remote_mapping
end
