function subgradient!(dd_model::DDModel)
    #get an upper bound
    alpha = dd_model.solver_data.alpha
    lower_bound = dd_model.current_lower_bound

    if dd_model.upper_bound_function != nothing
        dd_model.upper_bound = dd_model.upper_bound_function(dd_model)
        step = alpha*abs(lower_bound - dd_model.upper_bound)/(norm([dd_model.residual_equality;dd_model.residual_inequality])^2)
    else
        upper_bound = nothing
        step = alpha
    end

    lambda_equality_delta = step*dd_model.residuals_equality  #update multipliers
    lambda_inequality_delta = step*dd_model.residuals_inequality

    return lambda_equality_delta,lambda_inequality_delta
end

function cutting_plane!(dd_model::DDModel)
end
