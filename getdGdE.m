%% dGdEx, dGdEy, dGdEz Calculation
% This function calculates the derivatives of the objective function
% towards the components of the electric field at each grid point
function [dGdEx, dGdEy, dGdEz, dGda] = getdGdE(xv_dual, accelMatrix, ...
    Nx, Ny, Nz, ix_x, ix_y, ix_z, accelInterpMatrix)


    dGda = xv_dual' * (-accelMatrix);
    
    dGdax = dGda(ix_x);
    dGdEx = reshape(dGdax * accelInterpMatrix, [Nx Ny Nz]);

    dGday = dGda(ix_y);
    dGdEy = reshape(dGday * accelInterpMatrix, [Nx Ny Nz]);

    dGdaz = dGda(ix_z);
    dGdEz = reshape(dGdaz * accelInterpMatrix, [Nx Ny Nz]);
        
        
end