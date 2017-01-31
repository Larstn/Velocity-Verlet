%% Calculates the relative 2-norm of two values 
function [rel_error] = relative_n_norm(x1, x2,n)

    rel_error = abs((x1-x2)./x2 ).^n ;
    
end
