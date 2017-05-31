function [G_sum, dGdEx_sum, dGdEy_sum, dGdEz_sum, DG_sum, xv_dual] =...
    calculateDualVV(xv_all, objective_function, nParticle,Nt, ts, E_x, E_y, E_z, x_grid, y_grid, z_grid, n_charges, n_masses)

    elementary_charge   = 1.60217662e-19;
    electron_mass       = 1.6605e-27;
    Nx = size(E_x, 1);
    Ny = size(E_x, 2); 
    Nz = size(E_y, 3);
    
        
    G_sum = 0;
    dGdEx_sum = 0*E_x;
    dGdEy_sum = 0*E_y;
    dGdEz_sum = 0*E_z;
    DG_sum = zeros(6*Nt,1);
    
    [ix_x, ix_y, ix_z, ~] =...
        get_Index3D(Nt);
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

        [G, DG] = objective_function(xv);

        [systemMatrix, ~, accelMatrix] = ...
            velocityVerletMatrices3D(ts);

        systemMatrix = ...
            systemMatrix*...
            (1/((n_charges*elementary_charge)/(n_masses*electron_mass)));

        [S_p] = get_PrimalS(...
            Ix, Iy, Iz, E_x, E_y, E_z, Nt, ts, systemMatrix);

        xv_dual = S_p' \ DG;

        [dGdEx, dGdEy, dGdEz, ~] = getdGdE(xv_dual, accelMatrix, Nx, Ny, ...
            Nz, ix_x, ix_y, ix_z, accelInterpMatrix);
        
        G_sum = G_sum + G;
        
        DG_sum = DG_sum + 1/nParticle*2*DG*G;
        dGdEx_sum = dGdEx_sum + 1/nParticle*dGdEx;
        dGdEy_sum = dGdEy_sum + 1/nParticle*dGdEy;
        dGdEz_sum = dGdEz_sum + 1/nParticle*dGdEz;
        
    end
    
    

        G_sum = sqrt(1/nParticle*G_sum);
        dGdEx_sum = 0.5*dGdEx_sum ./ G_sum;
        dGdEy_sum = 0.5*dGdEy_sum ./ G_sum;
        dGdEz_sum = 0.5*dGdEz_sum ./ G_sum;

    end
        
        
        