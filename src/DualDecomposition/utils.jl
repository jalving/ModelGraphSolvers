
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

function prepare_link_matrix(link_constraints::Vector{LinkConstraint},sub_var_map::Dict)
    link_map = Dict()
    link_vars = []   #model node variables
    n_link_vars = 0
    I = Int64[]
    J = Int64[]
    V = Float64[]
    b = Float64[]
    i = 0
    for link in link_constraints  #Just enumerate here
        i += 1
        _unpack_constraint!(link,i,I,J,V,b,link_vars,link_map,sub_var_map)
    end
    A = sparse(I,J,V)

    vars = [sub_var_map[var] for var in link_vars]  #convert to subproblem variables
    return A,b,vars
end

function prepare_link_matrix(link_variables::Vector{LinkVariableRef},master_var_map::Dict)
    link_map = Dict()
    link_vars = []
    n_link_vars = 0
    I = Int64[]
    J = Int64[]
    V = Float64[]
    b = Float64[]
    row = 1
    for link_var in link_variables
        row += _unpack_constraint!(link_var,row,I,J,V,link_vars,link_map,master_var_map)
    end
    H = sparse(I,J,V)
    return H,link_vars #b is always zero for this matrix
end

function _unpack_constraint!(linkvariable::LinkVariableRef,row::Integer,I::Vector,J::Vector,V::Vector,link_vars::Vector,link_map::Dict,master_var_map::Dict)

    child_vars = master_var_map[linkvariable]
    n_sub_problems = length(child_vars)

    j = length(link_vars)
    for var in child_vars
        push!(link_vars,var)    #NOTE: There shouldn't be duplicate child vars from different link variables.      #A linked_variable should only have one parent
        j = length(link_vars)
        link_map[var] = j       #link index of variable var
    end

    for i = 1:n_sub_problems - 1

        var1 = child_vars[i]
        j1 = link_map[var1]             #look up the variable index since it was already added to the matrix

        var2 = child_vars[i+1]
        j2 = link_map[var2]             #look up the variable index since it was already added to the matrix

        coeff1 = -1
        coeff2 = 1

        push!(I,row)
        push!(J,j1)
        push!(V,coeff1)

        push!(I,row)
        push!(J,j2)
        push!(V,coeff2)

        row += 1
    end

    return row - 1
end


function _unpack_constraint!(constraint::LinkConstraint{JuMP.GenericAffExpr{Float64,JuMP.VariableRef},MOI.EqualTo{Float64}},row::Integer,I::Vector,J::Vector,V::Vector,b::Vector,link_vars::Vector,link_map::Dict,var_map::Dict)
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



#child_vars = linkvariable.graph.linkvariable_map[linkvariable]   #the children variables for this linkvariable

# j = length(link_vars)
# for var in child_vars
#     push!(link_vars,var)    #NOTE: There shouldn't be duplicate child vars from different link variables.      #A linked_variable should only have one parent
#     j = length(link_vars)
#     link_map[var] = j       #link index of variable var
# end

# for i = 1:length(child_vars) - 1    #for each child variable
#     var1 = child_vars[i]
#     j1 = link_map[var1]             #look up the variable index since it was already added to the matrix
#
#     var2 = child_vars[i+1]
#     j2 = link_map[var2]             #look up the variable index since it was already added to the matrix
#
#     coeff1 = -1
#     coeff2 = 1
#
#     push!(I,row)
#     push!(J,j1)
#     push!(V,coeff1)
#
#     push!(I,row)
#     push!(J,j2)
#     push!(V,coeff2)
#
#     row += 1
# end


#Add entry for master variable (link variable)
#master_var = linkvariable.vref                                  #the actual master JuMP variable
#push!(link_vars,master_var)
#coeff_master = 1
#j_master = length(link_vars)
#link_map[master_var] = j_master
