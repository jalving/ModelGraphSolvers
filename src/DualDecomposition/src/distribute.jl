#Distribute subproblems on a DD Model to workers.  return references to subproblems
function distribute_subproblems!(model::DDModel,to_workers::Vector{Int64})

    n_problems = length(model.subproblems)
    n_workers = length(to_workers)
    subproblems = model.subproblems
    problems_per_worker = Int64(floor(n_problems/n_workers))

    #Allocate modelnodes onto provided workers
    allocations = []
    j = 1
    while  j <= n_problems
        if j + problems_per_worker > n_problems
            push!(allocations,subproblems[j:end])
        else
            push!(allocations,subproblems[j:j+problems_per_worker - 1])
        end
        global j += problems_per_worker
    end

    worker_map = Dict(zip(workers,allocations))

    channel_map = Dict()
    for (worker,allocation) in worker_map
        channel = RemoteChannel(1)
        for subproblem in allocations
            @spawnat(1, put!(channel, subproblem))
        end
    end

    remote_subproblem_references = []
    remote_mapping = Dict()
    #Fill channel with subproblems
    @sync begin
        for
        # for (i,worker) in enumerate(to_workers)
            channel = worker_channel[worker]
            #@spawnat(1, put!(channel, allocations[i]))
            ref1 = @spawnat worker begin
                Core.eval(Main, Expr(:(=), :nodes, take!(channel)))
            end
            wait(ref1)
            push!(remote_references,ref1)
            remote_mapping[ref1] = worker
        end
        return remote_references
    end

    return remote_references, remote_mapping
end
