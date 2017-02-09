%% Input parameters for this test case

%Nt_vec = [2 8 20 40 100 200 400 800 1200 1600 2000 4000];
Nt_vec = [2 4 10 40 80 120];

Nx          = 21;
Ny          = 21;
Nz          = 2;

t_end = 1;
t_start = 0;

dGdEx_all = zeros(Nx, Ny, Nz, length(Nt_vec),2);

elementary_charge   = 1.60217662e-19;
electron_mass       = 9.1093856e-31;

xv_matrix = zeros(length(Nt_vec),6*Nt_vec(end));

dGdEx_diff = zeros(Nx,Ny,1,length(Nt_vec));

dGdEx_diff_fin = zeros(1,length(Nt_vec));
figure(200)
clf
figure(300)
clf

for i = 1:length(Nt_vec)


    Nt = Nt_vec(i);

    xv0 = [1.2; 1.2; 1.2; 0; 0; 0]*(elementary_charge/electron_mass);

    x_p  = 1.75*(elementary_charge/electron_mass);
    y_p  = 1.7*(elementary_charge/electron_mass);
    z_p  = 1.75*(elementary_charge/electron_mass);
    vx_p = 1*(elementary_charge/electron_mass);
    vy_p = 0*(elementary_charge/electron_mass);
    vz_p = 0*(elementary_charge/electron_mass);

    obj_weights = [1, 1, 1, 0, 0, 0];

    n_charges = 1;
    n_masses = 1;

    %% 
    % Creating test environment with given Input

    ts = linspace(t_start, t_end, Nt);

    x_grid = linspace(1, 2, Nx)*(elementary_charge/electron_mass);
    y_grid = linspace(1, 2, Ny)*(elementary_charge/electron_mass);
    z_grid = linspace(1, 2, Nz)*(elementary_charge/electron_mass);

    %V           = repmat(peaks(Nx), [1 1 2]); assert(Nx == Ny);

    %V(:,:,2) = V(:,:,1);

    E_x         = repmat(meshgrid(x_grid*0+1)', [1 1 2]);%-centeredDiff(V, 1);
    E_y         = repmat(meshgrid(y_grid*0+1),[1 1 2]);%-centeredDiff(V, 2);
    E_z         = repmat(meshgrid(0*x_grid), [1 1 2]);%-centeredDiff(V, 3);
    figure(1)
    subplot(3,1,1)
    imagesc(E_x(:,:,1))
    subplot(3,1,2)
    imagesc(E_y(:,:,1))
    subplot(3,1,3)
    imagesc(E_z(:,:,1))

    %%
    % Setting natural constants
    % in case that they are needed

    %elementary_charge   = 1.60217662e-19;
    %electron_mass       = 9.1093856e-31;
    %%


%     accelFunc = accelerationFunction( x_grid, y_grid, z_grid, ...
%         n_charges, n_masses);
    
    %[ix_x, ix_y, ix_z, iv_x, iv_y, iv_z, ia_x, ia_y, ia_z] = get_Index3D(Nt);
    
    

    nParticle = 1;

    hit_objective = @(xv) hitObjective3D_wrap(...
                xv, x_p, y_p, z_p, vx_p, vy_p, vz_p, obj_weights);

    [dGdEx0, dGdEy0, dGdEz0, G0, xv, DG, xv_dual, Nt] = ...
        VV_get_dual_E_final_auto_conv_4(n_charges, n_masses, E_x, E_y, E_z,...
    x_grid, y_grid, z_grid, xv0, nParticle, hit_objective, ts);



    %disp(G0)
    [ix_x, ix_y, ix_z, ~] =...
        get_Index3D(Nt);
    
    ts = linspace(t_start, t_end, Nt);
    x_end(i) = xv(max(ix_x));
   
     
   figure(10+i)
    subplot(3,15,1:2)
    plot(ts, xv(ix_x))
    title('x')

    subplot(3,15,16:17)
    plot(ts, xv(ix_y))
    title('y')
    
    
    subplot(3,15,31:32)
    plot(ts, xv(ix_z))
    title('z')
    
    xv_matrix(i,1:length(xv)) = xv;
    
    %%figure(10+i)
    subplot(3,15,[6:8 21:23 36:38] +1)
    imagesc(x_grid, y_grid, dGdEx0(:,:,1)');
    axis xy image
    colorbar
    title('dual')
    
    %figure(4)
    subplot(3,15,[9:11 24:26 39:40]+1)
    imagesc(x_grid, y_grid, dGdEx0(:,:,2)');
    axis xy image
    colorbar
    title('dual 2')
    
    dGdEx_all(:,:,:,i,2) = dGdEx0;
    
    
    figure(200)
    hold all 
    subplot(3,1,1)
    hold all
    plot(ts,xv(ix_x))
    subplot(3,1,2)
    hold all
    plot(ts,xv(ix_y))
    subplot(3,1,3)
    hold all
    plot(ts,xv(ix_z))
    
    accelFunc = accelerationFunction( x_grid, y_grid, z_grid, ...
    n_charges, n_masses);    
    [xv, accel] = velocityVerlet3D(ts, xv0, accelFunc(E_x, E_y, E_z));
    

    figure(300)
    hold all 
    subplot(3,1,1)
    hold all
    plot(ts,xv(ix_x))
    subplot(3,1,2)
    hold all
    plot(ts,xv(ix_y))
    subplot(3,1,3)
    hold all
    plot(ts,xv(ix_z))
    
    [G, dGdxv] = hitObjective3Dtarget(xv, obj_weights, x_p, y_p, z_p, ...
        vx_p, vy_p, vz_p);
    
    [systemMatrix, initMatrix, accelMatrix] = velocityVerletMatrices3D(ts);
    
    ax1 = accel(:,1);
    ay1 = accel(:,2);
    az1 = accel(:,3);
    
    delta = 1e-9;
    dGdEx_meas = 0*E_x(:,:,1);
%     
%     
%     for xx = 1:21
%         for yy = 1:21
%             fprintf('%i, %i\n', xx, yy);
% 
%             Ex2 = E_x;
%             Ex2(xx,yy,1) = Ex2(xx,yy,1) + delta;
%             [xv2, accel2] = velocityVerlet3D(ts, xv0, ...
%                 accelFunc(Ex2, E_y, E_z));
%             G2 = hitObjective3D(xv2, obj_weights);
%             G2b = hit_objective(xv2);
% 
%             disp(G2)
%             disp(G2b)
%             dGdEx_meas(xx,yy) = (G2b-G)/delta;
% 
% 
% 
% 
%         end
%     end
%     
        
    figure(10+i)
    subplot(3,15,[3:5 18:20 33:35] +1)
    imagesc(x_grid, y_grid, dGdEx_meas');
    axis xy image
    colorbar
    title('Meas')
    
    dGdEx_all(:,:,1,i,1) = dGdEx_meas;

    
    
    
    
    figure(20)
    subplot(length(Nt_vec),1,i)
    imagesc(x_grid, y_grid, (dGdEx0(:,:,1)' - dGdEx_meas')./dGdEx0(:,:,1)');
    colorbar
    
    
    figure(80)
    subplot(length(Nt_vec),1,i)
    imagesc(x_grid, y_grid, dGdEx0(:,:,1)');
    colorbar    
    
    
    dGdEx_diff(:,:,:,i) = (dGdEx0(:,:,1)' - dGdEx_meas')./dGdEx0(:,:,1)';
    indicies = dGdEx_diff(:,:,:,i);
    indicies(isnan(indicies)) = 0;
    dGdEx_diff_fin(i) = sum(sum(indicies));
    
    
    
   
    
end


% figure(3)
% plot(Nt_vec, dGdEx_diff_fin)



x_plot = zeros(1, length(Nt_vec)); 

for i = 1:length(Nt_vec)
    
    x_plot(i) = abs(x_end(end) - x_end(i)); 
    
end 

figure(6)
loglog(Nt_vec, x_plot)
xlabel('Nt_vec')
ylabel('x_{end_final} -  x_{end_i}')
title('x_end Convergance')


fin = zeros(1,length(Nt_vec));
for i = 1:length(Nt_vec)
    
    final_diff = (dGdEx_all(:,:,1,i,2)' - dGdEx_all(:,:,1,end,2)')./dGdEx_all(:,:,1,end,2)';
    final_diff(isnan(final_diff)) = 0;
    
    fin(i) = sum(sum(sqrt(final_diff.^2)));
    
end 
figure(4)
plot(Nt_vec,fin)
figure(5)
loglog(Nt_vec,fin)
xlabel('Nt_vec')
ylabel('Sum of Dual_{end} -  Sum of Dual_{end_i}')
title('Sum of Dual elements Convergance')
figure(3)
plot(Nt_vec, dGdEx_diff_fin)