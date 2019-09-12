
function get_link_eq_constraints(graph::ModelGraph)
    link_cons = LinkConstraint[]
    for link_con in getlinkconstraints(graph)
        if isa(link_con.set,MOI.EqualTo)
            push!(link_cons,link_con)
        end
    end
    return link_cons
end

function get_link_ineq_constraints(graph::ModelGraph)
    link_cons = LinkConstraint[]
    for link_con in getlinkconstraints(graph)
        if !(isa(link_con.set,MOI.EqualTo))
            push!(link_cons,link_con)  #It must be an inequality constraint
        end
    end
    return link_cons
end

function prepare_link_matrix(link_constraints::Vector{LinkConstraint})
    link_map = Dict()
    link_vars = []
    n_link_vars = 0
    I = Int64[]
    J = Int64[]
    V = Float64[]
    b = Float64[]
    i = 0
    for link in link_constraints  #Just enumerate here
        i += 1
        _unpack_constraint!(link,i,I,J,V,b,link_vars,link_map)
    end
    A = sparse(I,J,V)
    return A,b,link_vars,link_map
end

function prepare_link_matrix(link_variables::Vector{LinkVariableRef})
    link_map = Dict()
    link_vars = []
    n_link_vars = 0
    I = Int64[]
    J = Int64[]
    V = Float64[]
    b = Float64[]
    i = 0
    for link_var in link_variables  #Just enumerate here
        i += 1 #row
        _unpack_constraint!(link_var,i,I,J,V,link_vars,link_map)
    end
    A = sparse(I,J,V)
    return A,link_vars,link_map #b is always zero
end

function _unpack_constraint!(linkvariable::LinkVariableRef,row::Integer,I::Vector,J::Vector,V::Vector,link_vars::Vector,link_map::Dict)
    #Add entry for master variable (link variable)
    master_var = linkvariable.vref                                  #the actual master JuMP variable
    push!(link_vars,master_var)
    coeff = 1
    j = length(link_vars)
    push!(I,row)
    push!(J,j)
    push!(V,coeff)

    child_vars = linkvariable.graph.linkvariablemap[linkvariable]   #the children variables
    for var in child_vars
        if !(var in link_vars)
            j = length(link_vars) + 1
            link_map[var] = j       #link index of variable var
            push!(link_vars,var)
        else
            j = link_map[var]       #look up the variable index since it was already added to the matrix
        end
        coeff = -1
        push!(I,row)
        push!(J,j)
        push!(V,coeff)
    end
end


function _unpack_constraint!(constraint::LinkConstraint{JuMP.GenericAffExpr{Float64,JuMP.VariableRef},MOI.EqualTo{Float64}},row::Integer,I::Vector,J::Vector,V::Vector,b::Vector,link_vars::Vector,link_map::Dict)
    push!(b,constraint.set.value)
    for (var,coeff) in constraint.func.terms
        #check for new link constraint variables
        if !(var in link_vars)
            j = length(link_vars) + 1
            link_map[var] = j   #link index of variable var
            push!(link_vars,var)
        else
            #look up the variable index since it was already added to the matrix
            j = link_map[var]
        end
        push!(I,row)
        push!(J,j)
        push!(V,coeff)
    end
end

#NOTE: Might make this into a macro
function _unpack_constraint!(constraint::LinkConstraint{JuMP.GenericAffExpr{Float64,JuMP.VariableRef},MOI.LessThan{Float64}},row::Integer,I::Vector,J::Vector,V::Vector,b::Vector,link_vars::Vector,link_map::Dict)
    push!(b,constraint.set.upper)
    for (var,coeff) in constraint.func.terms
        #find a new link variable
        if !(var in link_vars)
            j = length(link_vars) + 1
            link_map[var] = j   #link index of variable var
            push!(link_vars,var)
        else
            #get correct variable index
            j = link_map[var]
        end
        push!(I,row)
        push!(J,j)
        push!(V,coeff)
    end
end

function _unpack_constraint!(constraint::LinkConstraint{JuMP.GenericAffExpr{Float64,JuMP.VariableRef},MOI.GreaterThan{Float64}},row::Integer,I::Vector,J::Vector,V::Vector,b::Vector,link_vars::Vector,link_map::Dict)
    push!(b,-1*constraint.set.lower)
    for (var,coeff) in constraint.func.terms
        #find a new link variable
        if !(var in link_vars)
            j = length(link_vars) + 1
            link_map[var] = j   #link index of variable var
            push!(link_vars,var)
        else
            #get correct variable index
            j = link_map[var]
        end
        push!(I,row)
        push!(J,j)
        push!(V,-1*coeff)
    end
end

#TODO: Interval Constraints
