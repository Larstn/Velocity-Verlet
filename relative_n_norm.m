%% Calculates the relative n-norm of two values 
function [rel_error] = relative_n_norm(x1, x2,n)

    rel_error = abs((x1-x2)./x2 ).^n ;
    rel_error(isnan(rel_error)) = 0;
end
