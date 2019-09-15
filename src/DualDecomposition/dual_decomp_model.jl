mutable struct DDSolverData
    max_iterations::Int64
    epsilon::Float64
    time_limit::Float64
    alpha::Float64
    delta::Float64
    max_no_improvement::Int64

    function DDSolverData()
        data = new()
        data.max_iterations = 100
        data.epsilon = 0.001
        data.time_limit = 3600
        data.alpha = 2.0
        data.delta = 0.5
        data.max_no_improvement = 10
        return data
    end
end

#IDEA: Create a DDModel using ModelGraph information
mutable struct DDModel
    solver_data::DDSolverData   #solver inputs
    solution::DDSolution        #solution data

    #Callbacks we can expose
    multiplier_update_function::Function                        #function that takes a DDModel argument and returns multiplier updates for equality and inequality multipliers
    upper_bound_function::Union{Nothing,Function}               #optional, can speed up convergence
    subproblem_solver::Union{Nothing,JuMP.OptimizerFactory}     #e.g. GLPK.Optimizer

    subproblems::Vector{JuMP.Model}                             #we assume these are all minimization problems

    linkvar_multipliers::Vector{Float64}
    equality_multipliers::Vector{Float64}                       #link equality constraint multipliers
    inequality_multipliers::Vector{Float64}                     #link inequality constraint multipliers

    link_var_matrix::Union{Nothing,AbstractMatrix}              #nonanticpitivity constraint matrix
    link_eq_matrix::Union{Nothing,AbstractMatrix}               #link equality matrix
    link_ineq_matrix::Union{Nothing,AbstractMatrix}             #link inequality matrix

    b_eq::Vector{Float64}                                       #link equality right hand side
    b_ineq::Vector{Float64}                                     #link inequality right hand side

    link_var_variables::Vector{JuMP.VariableRef}                #Vector of JuMP variables corresponding to nonanticpitivity constraints
    link_eq_variables::Vector{JuMP.VariableRef}                 #Vector of JuMP variables corresponding to link equality constraints
    link_ineq_variables::Vector{JuMP.VariableRef}               #Vector of JuMP variables corresponding to link inequality constraints

    current_lower_bound::Float64
    current_upper_bound::Float64

    residuals_linkvars::Vector           #A_linkvar*x_linkvar
    residuals_equality::Vector           #A_eq*x_eq - b_eq
    residuals_inequality::Vector         #A_ineq*x_ineq - b_ineq
end

# NOTE: We can lift shared variables into equality constraints
# shared_variables::Dict{JuMP.AbstractJuMPScalar,JuMP.AbstractJuMPScalar} #these get dualized into simple linkconstraints
# shared_variable_multipliers::Vector{Float64}
#dd_model.shared_variables = getsharedvariables(graph)  #dictionary

DDModel() = DDModel(DDSolverData(),DDSolution(),subgradient!,nothing,nothing,
                                Vector{JuMP.Model}(),
                                Vector{Float64}(),Vector{Float64}(),Vector{Float64}(),nothing,nothing,nothing,Vector{Float64}(),Vector{Float64}(),
                                Vector{JuMP.VariableRef}(),Vector{JuMP.VariableRef}(),Vector{JuMP.VariableRef}(),-Inf,Inf,
                                Vector{Float64}(),Vector{Float64}(),Vector{Float64}())

#Dual decomposition algorithm
#NOTE: This might get renamed to lagrange_solver depending how we deal with ADMM
function dual_decomposition_solve(graph::ModelGraph,args...;kwargs...)

    #NOTE: figure out which args to pass
    dd_model = DDModel(graph)
    status = solve(dd_model)

    #set solution on graph

    return status
end

function dual_decomposition_pmap_solve(graph::ModelGraph,workers,args...;kwargs...)
end

#Populate a DDModel using a ModelGraph
function DDModel(graph::ModelGraph) #args,kwargs

    #DDMODEL SUPPORTS
    #supports: linkconstraints
    #supports: linkvariables
    #does not support: subgraphs (need to aggregate)
    #does not support: interval constraints
    #does not support: nonlinear constraints

    #TODO Automatically aggregate subgraphs
    if has_subgraphs(graph)
        error("Dual Decomposition Solver does not support subgraphs.  run aggregate_subgraphs(graph) first" )
    end

    #Fill in all the data we need to do dual decomposition
    dd_model = DDModel()

    #Setup the subproblems based on modelnodes and linkvariables
    subproblems,sub_var_map,master_var_map =  build_subproblems(graph)

    dd_model.subproblems = subproblems               #copies model nodes and adds  master constraints to subproblems

    link_eq_constraints = get_link_eq_constraints(graph)          #link equality constraints in graph
    link_ineq_constraints = get_link_ineq_constraints(graph)      #link inequality constraints in graph
    link_variables = [var for var in getlinkvariables(graph)]     #link variables TODO: Make ModelGraphs return concrete type

    #Setup data structures
    link_eq_matrix,b_eq,link_eq_variables = prepare_link_matrix(link_eq_constraints,sub_var_map)
    link_ineq_matrix,b_ineq,link_ineq_variables, = prepare_link_matrix(link_ineq_constraints,sub_var_map)
    link_var_matrix,link_var_variables = prepare_link_matrix(link_variables,master_var_map)                 #This is the big sparse H matrix for nonanticpitivity constraints

    dd_model.link_var_matrix = link_var_matrix
    dd_model.link_var_variables = link_var_variables

    dd_model.link_eq_matrix = link_eq_matrix
    dd_model.b_eq = b_eq
    dd_model.link_eq_variables = link_eq_variables

    dd_model.link_ineq_matrix = link_ineq_matrix
    dd_model.b_ineq = b_ineq
    dd_model.link_ineq_variables = link_ineq_variables

    #Setup multipliers
    dd_model.equality_multipliers = zeros(size(link_eq_matrix,1))
    dd_model.inequality_multipliers = zeros(size(link_ineq_matrix,1))
    dd_model.linkvar_multipliers = zeros(size(link_var_matrix,1))

    #TODO Multiplier warm-start

    #Setup subproblem multiplier vectors
    _prepare_eq_multiplier_map!(link_eq_matrix,link_eq_variables)
    _prepare_ineq_multiplier_map!(link_ineq_matrix,link_ineq_variables)
    _prepare_linkvar_multiplier_map!(link_var_matrix,link_var_variables)

    return dd_model
end

#Solve a DD Model using dual decomposition or ADMM
function optimize!(dd_model::DDModel)

    println("Optimizing with Dual Decomposition:")
    println("# of subproblems: $(length(dd_model.subproblems))")
    println("# of nonanticpitivity constraints: ")
    println("# of link equality constraints:")
    println("# of link inequality constraints:")


    Zk_old = dd_model.current_lower_bound
    max_iterations = dd_model.solver_data.max_iterations
    epsilon = dd_model.solver_data.epsilon
    delta = dd_model.solver_data.delta

    no_improvement = 0
    for iter in 1:max_iterations
        println("Iteration: ",iter)
        # Solve subproblems
        Zk = 0  #objective value
        for subproblem in dd_model.subproblems
           Zkn = solve_subproblem!(dd_model,subproblem)
           Zk += Zkn                       # add objective values
        end
        #steptaken = true

        dd_model.current_lower_bound = Zk  #Update lower bound
        # If no improvement in the lowerbound, increase the no improvement counter
        if Zk < Zk_old
            no_improvement += 1
        end
        Zk_old = Zk

        # If too many iterations have happened without improvement in the lower bound, decrease alpha (step-size)
        if no_improvement >= dd_model.solver_data.max_no_improvement
            no_improvement = 0
            dd_model.solver_data.alpha *= delta
        end

        A_eq = dd_model.link_eq_matrix
        A_ineq = dd_model.link_ineq_matrix
        H_linkvars = dd_model.link_var_matrix

        x_eq = JuMP.value.(dd_model.link_eq_variables)
        x_ineq = JuMP.value.(dd_model.link_ineq_variables)
        x_linkvars = JuMP.value.(dd_model.link_var_variables)

        # Update residuals for multplier calculation
        dd_model.residuals_equality = A_eq*x_eq - dd_model.b_eq         #Ax = b
        dd_model.residuals_inequality = A_ineq*x_ineq - dd_model.b_ineq #Ax <= b
        dd_model.residuals_linkvars = H_linkvars*x_linkvars             #Hx = 0

        # Check convergence
        if norm(dd_model.residuals_equality) + norm(dd_model.residuals_linkvars) < epsilon && all(dd_model.residuals_inequality .>= 0)
            dd_model.solution.termination = :Optimal
            return :Optimal
        end

        #multiplier updates
        eq_multipliers_delta, ineq_multipliers_delta, linkvar_multipliers_delta = dd_model.multiplier_update_function(dd_model)
        dd_model.linkvar_multipliers += linkvar_multipliers_delta
        dd_model.equality_multipliers += eq_multipliers_delta
        dd_model.inequality_multipliers += ineq_multipliers_delta
        dd_model.inequality_multipliers[dd_model.inequality_multipliers .<= 0] .= 0  #NOTE: set to zero, or just don't update?
    end

    dd_model.solution.termination= :MaxIterationsReached
    return :MaxIterationsReached

end

#Use current multiplier values to solve a lagrange sub-problem
function solve_subproblem!(dd_model::DDModel,subproblem::JuMP.Model)

    #IDEA: I could use parameter JuMP and make the multipliers into parameters
    m = subproblem

    eq_multipliers = dd_model.equality_multipliers
    ineq_multipliers = dd_model.inequality_multipliers
    linkvar_multipliers = dd_model.linkvar_multipliers

    scale = m.ext[:objective_scale]
    JuMP.set_objective_function(m,scale*m.ext[:original_objective])

    #Add dualized part to objective function for each linkconstraint that hits this node
    obj = JuMP.objective_function(m)

    #Add equality multiplier terms
    for term in m.ext[:eq_multipliers]  #could have duplicate terms
        coeff = term[1]
        var = term[2]
        multiplier_index = term[3]
        obj += coeff*var*eq_multipliers[multiplier_index]
    end

    #Add inequality multiplier terms
    for term in m.ext[:ineq_multipliers]
        coeff = term[1]
        var = term[2]
        multiplier_index = term[3]
        obj += coeff*var*ineq_multipliers[multiplier_index]
    end

    #Add nonanticpitivity multiplier terms
    for term in m.ext[:linkvar_multipliers]
        coeff = term[1]
        var = term[2]
        multiplier_index = term[3]
        obj += coeff*var*linkvar_multipliers[multiplier_index]
    end

    JuMP.set_objective_function(m,obj)

    JuMP.optimize!(m)

    m.ext[:current_objective_value] = JuMP.objective_value(m)
    return JuMP.objective_value(m)
end


function _prepare_linkvar_multiplier_map!(link_matrix::SparseMatrixCSC{Float64,Int64},link_variables::Vector)
    for i in 1:size(link_matrix)[1] #loop through rows
        row = link_matrix[i,:]
        for j in row.nzind
            coeff = link_matrix[i,j]
            var = link_variables[j]
            model = var.model
            push!(model.ext[:linkvar_multipliers],(coeff,var,i))
        end
    end
end


function _prepare_eq_multiplier_map!(link_matrix::SparseMatrixCSC{Float64,Int64},link_variables::Vector)
    for i in 1:size(link_matrix)[1] #loop through rows
        row = link_matrix[i,:]
        for j in row.nzind
            coeff = link_matrix[i,j]
            var = link_variables[j]
            model = var.model
            push!(model.ext[:eq_multipliers],(coeff,var,i))
        end
    end
end

function _prepare_ineq_multiplier_map!(link_matrix::SparseMatrixCSC{Float64,Int64},link_variables::Vector)
    for i in 1:size(link_matrix)[1] #loop through rows
        row = link_matrix[i,:]
        for j in row.nzind
            coeff = link_matrix[i,j]
            var = link_variables[j]
            model = var.model
            push!(model.ext[:ineq_multipliers],(coeff,var,i))
        end
    end
end



function build_subproblems(graph::ModelGraph)

    linkvariables = getlinkvariables(graph)
    sub_var_map = Dict()
    master_var_map = Dict((k,[]) for k in linkvariables)


    subproblems = []
    for node in getnodes(graph)
        model = getmodel(node)
        model.ext = Dict()          #need to clear out the model.ext in order to copy

        subproblem = copy(model)
        #Attach optimizer from node

        ref_map = ModelGraphs.AggregationMap(subproblem)

        #ref_map.varmap =  Dict(zip(JuMP.all_variables(model),JuMP.all_variables(subproblem)))
        new_vars = JuMP.all_variables(subproblem)
        for var in JuMP.all_variables(model)
            index = var.index
            new_var = new_vars[index.value]
            ref_map[var] = new_vars[index.value]
            if var in keys(node.linkvariablemap)     #Point linked_variables to master problem variables
                push!(master_var_map[node.linkvariablemap[var]],new_vars[index.value])
                ref_map[node.linkvariablemap[var].vref] = new_vars[index.value]
            end
        end

        subproblem.ext[:modelnode] = node

        #Add extra link variables if necessary to submodels.
        for linkvariable in linkvariables
            if !(linkvariable in values(node.linkvariablemap)) #if the child does not have a master variable
                new_var = @variable(subproblem)
                var_name = JuMP.name(linkvariable.vref)
                new_name = var_name
                JuMP.set_name(new_var,new_name)
                if JuMP.start_value(linkvariable.vref) != nothing
                    JuMP.set_start_value(new_x,JuMP.start_value(linkvariable.vref))
                end
                # ref_map[new_var] = linkvariable.vref
                ref_map[linkvariable.vref] = new_var
            else

            end
        end


        constraint_types = JuMP.list_of_constraint_types(graph.mastermodel)
        for (func,set) in constraint_types
            constraint_refs = JuMP.all_constraints(graph.mastermodel, func, set)

            #Don't add duplicate constraints for master variable
            # if func == JuMP.VariableRef
            #     var = constraint.func
            #     if is_linked_variable(var)  #Don't add multiple copies of a linked variable constraint
            #         continue
            #     end
            # end


            for constraint_ref in constraint_refs
                constraint = JuMP.constraint_object(constraint_ref)
                new_constraint = ModelGraphs._copy_constraint(constraint,ref_map)
                new_ref = JuMP.add_constraint(subproblem,new_constraint)
                ref_map[constraint_ref] = new_ref
            end
        end

        subproblem.ext[:eq_multipliers] = []
        subproblem.ext[:ineq_multipliers] = []
        subproblem.ext[:linkvar_multipliers] = []
        subproblem.ext[:objective_scale] = 1
        obj = JuMP.objective_function(subproblem)
        if JuMP.objective_sense(subproblem) == MOI.MAX_SENSE
            subproblem.ext[:objective_scale] = -1
            JuMP.set_objective_sense(subproblem,MOI.MIN_SENSE)
        end

        subproblem.ext[:original_objective] = obj
        push!(subproblems,subproblem)

        merge!(sub_var_map,ref_map.varmap)


    end
    return subproblems, sub_var_map, master_var_map
end

# #Flip objective, add dictionary entry for multipliers for this model
# function prepare_subproblem!(graph::ModelGraph,m::JuMP.Model)
#     m.ext[:eq_multipliers] = []
#     m.ext[:ineq_multipliers] = []
#     m.ext[:objective_scale] = 1
#     obj = JuMP.objective_function(m)
#     if JuMP.objective_sense(m) == MOI.MAX_SENSE
#         m.ext[:objective_scale] = -1
#         JuMP.set_objective_sense(m,MOI.MIN_SENSE)
#     end
#     m.ext[:original_objective] = obj
#     #Add master constraints corresponding to linked_variable
# end
