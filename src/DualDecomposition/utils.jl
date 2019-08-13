
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
    #link_constraints = graph.linkmodel.linkconstraints #index => constraint
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

function _unpack_constraint!(constraint::LinkConstraint{JuMP.GenericAffExpr{Float64,JuMP.VariableRef},MOI.EqualTo{Float64}},row::Integer,I::Vector,J::Vector,V::Vector,b::Vector,link_vars::Vector,link_map::Dict)
    push!(b,constraint.set.value)
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
#TODO: Interval

# function prepare_link_ineq_matrix(graph::ModelGraph)
#     #link_constraints = graph.linkmodel.linkconstraints #index => constraint
#     link_map = Dict()
#     link_vars = []
#     n_link_vars = 0
#     I = Int64[]
#     J = Int64[]
#     V = Float64[]
#     b = Float64[]
#     i = 0
#     for link_con in getlinkconstraints(graph)
#         #If it's an equality constraint
#
#         if isa(link_con.set,MOI.LessThan)  #Good to go
#             i += 1
#
#             push!(b,link_con.set.upper)
#             #Don't need to modify coefficients
#             for (var,coeff) in link_con.func.terms
#                 #find a new link variable
#                 if !(var in keys(link_vars))
#                     n_link_vars += 1
#                     j = n_link_vars
#                     link_map[var] = j   #link index of variable var
#                     push!(link_vars,var)
#                 else
#                     #get correct variable index
#                     j = link_vars[var]
#                 end
#
#                 push!(I,i)
#                 push!(J,j)
#                 push!(V,coeff)
#             end
#
#         elseif isa(link_con.set,MOI.GreaterThan)
#             i += 1
#             push!(b,-1*link_con.set.lower)
#             #Don't need to modify coefficients
#             for (var,coeff) in link_con.func.terms
#                 #find a new link variable
#                 if !(var in keys(link_vars))
#                     n_link_vars += 1
#                     j = n_link_vars
#                     link_map[var] = j   #link index of variable var
#                     push!(link_vars,var)
#                 else
#                     #get correct variable index
#                     j = link_vars[var]
#                 end
#
#                 push!(I,i)
#                 push!(J,j)
#                 push!(V,-1*coeff)
#             end
#
#         elseif isa(link_con.set,MOI.Interval)
#             error("Lagrange Solver does not yet support Interval constraints")
#         end
#     end
#
#     Pi = sparse(I,J,V)
#     return Pi,link_vars,b
# end
