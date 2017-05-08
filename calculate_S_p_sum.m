%% Calculate S_primal for a sum objective Function, e.g. RMS

function [S_p_sum, accelInterpMatrix_sum, accelMatrix] = calculate_S_p_sum(xv_all, x_grid, y_grid, z_grid,...
   Nx, Ny, Nz, ts, n_charges, elementary_charge, n_masses, electron_mass, E_x, E_y, E_z) 

    Nt = size(xv_all,1)./6;
    nParticle = size(xv_all,2);
    
    [ix_x, ix_y, ix_z, ~] =...
                     get_Index3D(Nt);
    
     S_p_sum = zeros(6*Nt, 6*Nt);
     accelInterpMatrix_sum = zeros(Nt, Nx*Ny*Nz);
    
    for i = 1:nParticle
        
        xv = xv_all(:,i);
        [iix, iiy, iiz, w000, w001, w010, w011, w100, w101, w110, w111] ...
                = trilinear_weights(xv(ix_x), xv(ix_y), xv(ix_z), ...
                    x_grid, y_grid, z_grid);


        accelInterpMatrix = get_accelInterpmatrix(iix, iiy, iiz, ...
            w000, w001, w010, w011, w100, w101, w110, w111, ...
            Nx, Ny, Nz, Nt);

        [Ix, Iy, Iz] = get_I(xv, ix_x, ix_y, ix_z,...
            x_grid, y_grid, z_grid, Nx, Ny, Nz, Nt);


        [systemMatrix, ~, accelMatrix] = ...
            velocityVerletMatrices3D(ts);

        systemMatrix = ...
            systemMatrix*...
            (1/((n_charges*elementary_charge)/(n_masses*electron_mass)));

        [S_p] = get_PrimalS(...
            Ix, Iy, Iz, E_x, E_y, E_z, Nt, ts, systemMatrix);


        S_p_sum = S_p_sum + S_p;
        accelInterpMatrix_sum = accelInterpMatrix_sum + accelInterpMatrix;
    end
end 


