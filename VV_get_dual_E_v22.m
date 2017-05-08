%% Calculates dGdE(x,y,z)
% Calculates the derivative of the objective function with respect to the
% electric field components at each grid with the dual function
function [dGdEx_sum, dGdEy_sum, dGdEz_sum, G_sum, xv_all, DG_sum, xv_dual, Nt] ...
    = VV_get_dual_E_v22(n_charges, n_masses, V_Comsol, ...
    x_grid, y_grid, z_grid, xv0, nParticle, objective_function, ts)

    E_x = V_Comsol(:,:,:,1);
    E_y = V_Comsol(:,:,:,2);
    E_z = V_Comsol(:,:,:,3);
    
    elementary_charge   = 1.60217662e-19;
    electron_mass       = 1.6605e-27;

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

     
    for ii = 1:nParticle
        
        cnt = 0;
        cnt_no_brakes = 0;
        xv = zeros(6*Nt, 1);

        while(1)

            xv_start = xv0(:,ii);
            xv_prev = xv;
            ix_x_prev = ix_x;
            ix_y_prev = ix_y;
            ix_z_prev = ix_z;
            
            [ix_x, ix_y, ix_z, ~] =...
                     get_Index3D(Nt);

            [xv, ~] = velocityVerlet3D(...
                ts, xv_start, accelFunc(E_x, E_y, E_z));
 
              if cnt > 0
                diff_x = ...
                    relative_n_norm(xv(ix_x(end)),xv_prev(ix_x_prev(end)),2)
                diff_y = ...
                    relative_n_norm(xv(ix_y(end)),xv_prev(ix_y_prev(end)),2)
                diff_z = ...
                    relative_n_norm(xv(ix_z(end)), xv_prev(ix_z_prev(end)),2)

                
                delta_x = max(abs(xv(ix_x(1:(end-1))) - xv(ix_x(2:end))));
                delta_y = max(abs(xv(ix_y(1:(end-1))) - xv(ix_y(2:end)))); 
                delta_z = max(abs(xv(ix_z(1:(end-1))) - xv(ix_z(2:end))));
                
                

                if (diff_x < 1e-14) && (diff_y < 1e-14) &&...
                        (diff_z < 1e-14) && (delta_x < 0.5*d_x_diff) ...
                        && (delta_y < 0.5*d_y_diff) && (delta_z < 0.5*d_z_diff)
                    disp('all good')
                    Nt = Nt/2;
                    ts = linspace(t_start,t_end,Nt);
                    break                  
                else
                    disp('not good')
                    Nt = 2*Nt;
                    ts = linspace(t_start,t_end,Nt);
                    cnt = cnt + 1;
                    
                end
                
              elseif cnt == 0
                  
                    disp('f?rste gang')
                    Nt = 2*Nt;
                    ts = linspace(t_start,t_end,Nt);
                    cnt = cnt + 1;
                    continue
               end
                   
        end
      
    

        
        
    end
    xv_all = zeros(6*Nt,nParticle);
    [ix_x, ix_y, ix_z, ~] =...
                     get_Index3D(Nt);
                 
  disp('Final number of time points:')
  disp(Nt)

  for i = 1:nParticle
      [xv, ~] = velocityVerlet3D(...
        ts, xv_start, accelFunc(E_x, E_y, E_z));
        xv_all(:,i) = xv;
  end 

        
                [S_p_sum, accelInterpMatrix_sum, accelMatrix] = calculate_S_p_sum(xv_all, x_grid, y_grid, z_grid,...
   Nx, Ny, Nz, ts, n_charges, elementary_charge, n_masses, electron_mass, E_x, E_y, E_z);
    

        
        [G_sum, DG_sum] = objective_function(xv_all);

        xv_dual = S_p_sum' \ DG_sum;

        [dGdEx_sum, dGdEy_sum, dGdEz_sum, ~] = getdGdE(xv_dual, accelMatrix, Nx, Ny, ...
            Nz, ix_x, ix_y, ix_z, accelInterpMatrix_sum);
        

end

