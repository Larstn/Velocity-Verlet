%% Calculates dGdE(x,y,z)
% Calculates the derivative of the objective function with respect to the
% electric field components at each grid with the dual function
function [dGdEx_sum, dGdEy_sum, dGdEz_sum, G_sum, xv_all, DG_sum, xv_dual, Nt] ...
    = VV_get_dual_E_v21(n_charges, n_masses, E_x, E_y, E_z, ...
    x_grid, y_grid, z_grid, xv0, nParticle, objective_function, ts)


    elementary_charge   = 1.60217662e-19;
    electron_mass       = 1.6605e-27;

    Nx = size(E_x, 1);
    Ny = size(E_x, 2); 
    Nz = size(E_y, 3);
    
    accelFunc = accelerationFunction( x_grid, y_grid, z_grid, ...
    n_charges, n_masses);

    %Nt = 100;
   % ts = linspace(t_start, t_end, Nt); % second
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
    
%     disp('d_x_diff')
%     disp(d_x_diff)
%     disp('d_y_diff')
%     disp(d_y_diff)
%     disp('d_z_diff')
%     disp(d_z_diff)
    xv   = zeros(6*Nt, nParticle);
    
    
     G_sum = 0;
%     G_sum_old = 0;
    dGdEx_sum = 0*E_x;
    dGdEy_sum = 0*E_y;
    dGdEz_sum = 0*E_z;
     DG_sum = zeros(6*Nt,1)

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
        %xv
        
        xv_all(:,ii) = xv;
    
    
    
     [iix, iiy, iiz, w000, w001, w010, w011, w100, w101, w110, w111] ...
            = trilinear_weights(xv(ix_x), xv(ix_y), xv(ix_z), ...
                x_grid, y_grid, z_grid);


        accelInterpMatrix = get_accelInterpmatrix(iix, iiy, iiz, ...
            w000, w001, w010, w011, w100, w101, w110, w111, ...
            Nx, Ny, Nz, Nt);

        [Ix, Iy, Iz] = get_I(xv, ix_x, ix_y, ix_z,...
            x_grid, y_grid, z_grid, Nx, Ny, Nz, Nt);

        [G, DG] = objective_function(xv);
        %G_sum     = G_sum + G.^2;
        %G_sum_old = G_sum_old + G;
       % DG_sum = DG_sum + 1/nParticle*2*DG*G;
        %G_sum = G_sum + G;
        %DG_sum = DG_sum + DG;

   
    %G_sum = sqrt(1/nParticle*G_sum);
    %DG_sum = 1/(2*G_sum) * DG_sum;
    %DG_sum = G_sum_old * 2 * DG_sum * 1/(sqrt(G_sum_old^2)*2);
    
    
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
       % DG = 0;
       % xv_dual = 0;
            
    end
    

    G_sum = sqrt(1/nParticle*G_sum);
    dGdEx_sum = 0.5*dGdEx_sum ./ G_sum;
    dGdEy_sum = 0.5*dGdEy_sum ./ G_sum;
    dGdEz_sum = 0.5*dGdEz_sum ./ G_sum;

end

