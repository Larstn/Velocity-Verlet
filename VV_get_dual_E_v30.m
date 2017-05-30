%% Calculates dGdE(x,y,z)
% Calculates the derivative of the objective function with respect to the
% electric field components at each grid with the dual function
function [dGdEx_sum, dGdEy_sum, dGdEz_sum, G_sum, xv_all, DG_sum, xv_dual, Nt, rms_r_spot] ...
    = VV_get_dual_E_v30(n_charges, n_masses, V_Comsol, ...
    x_grid, y_grid, z_grid, xv0, nParticle, objective_function, ts, k, x_p, y_p, z_p)

    E_x = V_Comsol(:,:,:,1); 
    E_y = V_Comsol(:,:,:,2);
    E_z = V_Comsol(:,:,:,3);

    


    Nx = size(E_x, 1);
    Ny = size(E_x, 2); 
    Nz = size(E_y, 3);
    
    accelFunc = accelerationFunction( x_grid, y_grid, z_grid, ...
    n_charges, n_masses);

    Nt = length(ts);
    t_start = ts(1);
    t_end = ts(end);
    
    [ix_x, ix_y, ix_z, ~] =...
        get_Index3D(Nt);

    
    d_x_diff = x_grid(2) - x_grid(1);
    d_y_diff = y_grid(2) - y_grid(1);
    if length(z_grid) == 1
        d_z_diff = 0;
    else
        d_z_diff = z_grid(2) - z_grid(1);
    end
    
    xv   = zeros(6*Nt, nParticle);
    


    for ii = 1:nParticle
        
        cnt = 0;
        
        while(1)
            cnt = cnt+1;
            xv_start = xv0(:,ii);
            xv_prev = xv;
            ix_x_prev = ix_x;
            ix_y_prev = ix_y;
            ix_z_prev = ix_z;
            
            [ix_x, ix_y, ix_z, ~] =...
                     get_Index3D(Nt);

            [xv, ~] = velocityVerlet3D(...
                ts, xv_start, accelFunc(E_x, E_y, E_z));
%                 if cnt == 1
%                     break
%                 end
               if cnt > 0 
                
                diff_x = ...
                    relative_n_norm(xv(ix_x(end)),xv_prev(ix_x_prev(end)),2);
                diff_y = ...
                    relative_n_norm(xv(ix_y(end)),xv_prev(ix_y_prev(end)),2);
                diff_z = ...
                    relative_n_norm(xv(ix_z(end)), xv_prev(ix_z_prev(end)),2);
                disp('Nt:')
                disp(Nt)
%                 disp('diff_x:')
%                 disp(diff_x)
%                 disp(diff_y)
%                 disp(diff_z)
                
                delta_x = max(abs(xv(ix_x(1:(end-1))) - xv(ix_x(2:end))));
                delta_y = max(abs(xv(ix_y(1:(end-1))) - xv(ix_y(2:end)))); 
                delta_z = max(abs(xv(ix_z(1:(end-1))) - xv(ix_z(2:end)))); 
                
%                 disp('delta x y z')
%                 disp(delta_x)
%                 disp(delta_y)
%                 disp(delta_z)

                if (diff_x < 0.0001) && (diff_y < 0.0001) &&...
                        (diff_z < 0.0001) && (delta_x < 0.5*d_x_diff) ...
                        && (delta_y < 0.5*d_y_diff) && (delta_z < 0.5*d_z_diff)
                    break
                elseif (ii >1) 
                    break                  
                elseif cnt == 1
                    break
                else
                    Nt = 2*Nt;
                    ts = linspace(t_start,t_end,Nt);
                    xv_all    = zeros(6*Nt, nParticle);
                    
                end
                
                
              
                
               end
                   
        end
        
        
        xv_all(:,ii) = xv;
    
    end
    
    fprintf('Done Trajectories \n');
    
    [G_sum, dGdEx_sum, dGdEy_sum, dGdEz_sum, DG_sum, xv_dual] = calculateDualVV(xv_all, ...
        objective_function, nParticle, Nt, ts, E_x, E_y, E_z,...
        x_grid, y_grid, z_grid, n_charges, n_masses);
    
    
    
    [z_spot_fin, y_spot_fin, rms_r_spot] = calculateSpotsize(xv_all, Nt, V_Comsol(:,:,:,1), V_Comsol(:,:,:,2), V_Comsol(:,:,:,3), x_grid, y_grid, z_grid, k,ts, n_charges, n_masses, nParticle, x_p, y_p, z_p);


end

