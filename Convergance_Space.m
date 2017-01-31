%% Input parameters for this test case
Nxy_vec = [2 5 10 20 40 80 160 320];

Nx          = 11;
Ny          = 11;
Nz          = 2;

t_max = 3;

dGdEx_all = zeros(Nx, Ny, Nz, length(Nt_vec));

elementary_charge   = 1.60217662e-19;
electron_mass       = 9.1093856e-31;

for i = 1:length(Nxy_vec)

    Nx = Nxy_vec(i);
    Ny = Nx;
    Nt = 100;

    xv0 = [0; 0; 0; .8; .499; 0]*(elementary_charge/electron_mass);

    x_p = 0.25*(elementary_charge/electron_mass);
    y_p  = 0.3*(elementary_charge/electron_mass);
    z_p  = 0*(elementary_charge/electron_mass);
    vx_p = 1*(elementary_charge/electron_mass);
    vy_p = 0*(elementary_charge/electron_mass);
    vz_p = 0*(elementary_charge/electron_mass);

    obj_weights = [1, 1, 1, 0, 0, 0];

    n_charges = 1;
    n_masses = 1;

    %% 
    % Creating test environment with given Input

    ts = linspace(0, 2, Nt);

    x_grid = linspace(-1, 1, Nx)*(elementary_charge/electron_mass);
    y_grid = linspace(-1, 1, Ny)*(elementary_charge/electron_mass);
    z_grid = linspace(-1, 1, Nz)*(elementary_charge/electron_mass);

    V           = repmat(peaks(Nx), [1 1 2]); assert(Nx == Ny);

    V(:,:,2) = V(:,:,1);

    E_x         = -centeredDiff(V, 1);
    E_y         = -centeredDiff(V, 2);
    E_z         = -centeredDiff(V, 3);
    figure(2)
    subplot(3,1,1)
    imagesc(E_x(:,:,1))
    subplot(3,1,2)
    imagesc(E_y(:,:,1))
    subplot(3,1,3)
    imagesc(E_z(:,:,1))

    %%
    % Setting natural constants
    % in case that they are needed

    elementary_charge   = 1.60217662e-19;
    electron_mass       = 9.1093856e-31;
    %%


    accelFunc = accelerationFunction( x_grid, y_grid, z_grid, ...
        n_charges, n_masses);
    
    %[ix_x, ix_y, ix_z, iv_x, iv_y, iv_z, ia_x, ia_y, ia_z] = get_Index3D(Nt);
    
    [xv, accel] = velocityVerlet3D(ts, xv0, accelFunc(E_x, E_y, E_z));
    % 
    % [G, dGdxv] = hitObjective3Dtarget(xv, obj_weights, x_p, y_p, z_p, ...
    %     vx_p, vy_p, vz_p);
    % 
    % [systemMatrix, initMatrix, accelMatrix] = velocityVerletMatrices3D(ts);
    % 
    % ax1 = accel(:,1);
    % ay1 = accel(:,2);
    % az1 = accel(:,3);
    % 
    % delta = 1e-9;
    % dGdEx_meas = 0*E_x(:,:,1);
    % 
    % 
    % for xx = 5:8
    % for yy = 5:8
    %     fprintf('%i, %i\n', xx, yy);
    %     
    %     Ex2 = E_x;
    %     Ex2(xx,yy,1) = Ex2(xx,yy,1) + delta;
    %     [xv2, accel2] = velocityVerlet3D(ts, xv0, accelFunc(Ex2, E_y, E_z));
    %     G2 = hitObjective3D(xv2, obj_weights);
    % 
    %     
    %     dGdEx_meas(xx,yy) = (G2-G)/delta;
    %     
    % 
    %     
    %     
    % end
    % end
    % 
    % figure(4)
    % subplot(2,1,2)
    % imagesc(x_grid, y_grid, dGdEx_meas');
    % axis xy image
    % colorbar
    % title('Meas')


    nParticle = 1;

    hit_objective = @(xv) hitObjective3D_wrap(...
                xv, x_p, y_p, z_p, vx_p, vy_p, vz_p, obj_weights);

    [dGdEx0, dGdEy0, dGdEz0, G0, xv, DG, xv_dual] = VV_get_dual_E_final(n_charges, n_masses, E_x, E_y, E_z,  x_grid, y_grid, z_grid, xv0, nParticle, hit_objective);

    [ix_x, ix_y, ix_z, ~] =...
        get_Index3D(Nt);
     ts = linspace(0, 2, Nt);
    figure(5)
    subplot(3,1,1)
    plot(ts, xv(ix_x))

    subplot(3,1,2)
    plot(ts, xv(ix_y))
    
    
    subplot(3,1,3)
    plot(ts, xv(ix_z))
    
    
    figure(4)
    subplot(length(Nt_vec),2,2*i - 1)
    imagesc(x_grid, y_grid, dGdEx0(:,:,1)');
    axis xy image
    colorbar
    title('dual')
    
    figure(4)
    subplot(length(Nt_vec),2,2*i)
    imagesc(x_grid, y_grid, dGdEx0(:,:,2)');
    axis xy image
    colorbar
    title('dual')


%     Nx          = 11;
%     Ny          = 11;
%     Nz          = 2;
% 
% 
% 
%     V           = repmat(peaks(Nx), [1 1 2]); assert(Nx == Ny);
%     V(:,:,2) = V(:,:,1);
%     E_x         = -centeredDiff(V, 1);
%     E_y         = -centeredDiff(V, 2);
%     E_z         = -centeredDiff(V, 3);
%     T_span      = [0 1];
% 
%     dGdEx_all(:,:,:,1) = dGdEx0;

    %sum(sum(abs(dGdEx0(:,:,1)' - dGdEx_meas')))
    %sum(abs(dGdEx0(:,:,1)' - dGdEx_meas).^2)
    
end