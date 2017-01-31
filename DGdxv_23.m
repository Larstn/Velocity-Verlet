%% Derivative towards the other directions
% Calculates the derivative of G towards the elements of the imput vector
% that are not defining the plane

function [DGdxv1, DGdxv2] = DGdxv_23(w ,xv_relevant, wobj_1, wobj_2, ...
    xv, xv_p)

    
    DGdxv1 = ...
        w*2*(wobj_1*xv(xv_relevant(1))...
        + wobj_2*xv(xv_relevant(2)) - xv_p)*wobj_1;
    
    DGdxv2 = ...
       w*2*(wobj_1*xv(xv_relevant(1))...
       + wobj_2*xv(xv_relevant(2)) - xv_p)*wobj_2;
    
    

end