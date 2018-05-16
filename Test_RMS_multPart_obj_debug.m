%% Input parameters for this test case
%clc 
%clear all
%%
% Setting natural constants
% in case that they are needed

elementary_charge   = 1.60217662e-19;
electron_mass       = 9.1093856e-31;

Nx          = 11;
Ny          = 11;
Nz          = 2;

Nt = 100;

x_grid = linspace(-1*(elementary_charge/electron_mass), 1*(elementary_charge/electron_mass), Nx);
y_grid = linspace(-1*(elementary_charge/electron_mass), 1*(elementary_charge/electron_mass), Ny);
z_grid = linspace(-1*(elementary_charge/electron_mass), 1*(elementary_charge/electron_mass), Nz);

xv0 = [0 0; 0 0; 0 0;.8 .8; .499 .499; 0 0]*(elementary_charge/electron_mass);

x_p = 0.25*(elementary_charge/electron_mass);
y_p  = 0.3*(elementary_charge/electron_mass);
z_p  = 0*(elementary_charge/electron_mass);
vx_p = 1*(elementary_charge/electron_mass);
vy_p = 0*(elementary_charge/electron_mass);
vz_p = 0*(elementary_charge/electron_mass);

obj_weights = [1, 1, 1, 0, 0, 0];

n_charges = 1;
n_masses = 1;
nParticle = 2;

%% 
% Creating test environment with given Input

ts = linspace(0, 1, Nt);



V           = repmat(peaks(Nx), [1 1 2]); assert(Nx == Ny);

V(:,:,2) = V(:,:,1);

E_x         = -centeredDiff(V, 1);
E_y         = -centeredDiff(V, 2);
E_z         = -centeredDiff(V, 3);


%%
    figure(20)
    subplot(3,1,1)
    imagesc(E_x(:,:,1))
    subplot(3,1,2)
    imagesc(E_y(:,:,1))
    subplot(3,1,3)
    imagesc(E_z(:,:,1))
    
%%    

hit_objective = @(x_v) hitObjective3D_wrap(...
            x_v, x_p, y_p, z_p, vx_p, vy_p, vz_p, obj_weights);


%%
dGdEx_meas = 0*E_x(:,:,1);

accelFunc = accelerationFunction( x_grid, y_grid, z_grid, ...
    n_charges, n_masses);


[xv, accel] = velocityVerlet3D(ts, xv0(:,1), accelFunc(E_x, E_y, E_z));
[xv2, accel] = velocityVerlet3D(ts, xv0(:,2), accelFunc(E_x, E_y, E_z));


G = sqrt(1/nParticle*(hit_objective(xv).^2 + hit_objective(xv2).^2));

delta = 1e-3;


for xx = 4:9%1:length(x_grid)
    for yy = 4:9%1:length(y_grid)
        fprintf('%i, %i\n', xx, yy);

        Ex2 = E_x;
        Ex2(xx,yy,1) = Ex2(xx,yy,1) + delta;
        [xv3, accel2] = velocityVerlet3D(ts, xv0, accelFunc(Ex2, E_y, E_z));
        [xv4, accel2] = velocityVerlet3D(ts, xv0, accelFunc(Ex2, E_y, E_z));

        G2 = sqrt(1/nParticle*(hit_objective(xv3).^2 + hit_objective(xv4).^2));

        dGdEx_meas(xx,yy) = dGdEx_meas(xx,yy) + (G2-G)/delta;




    end
end

% 
% dGdEx_meas2 = 0*E_x(:,:,1);
% 
% accelFunc = accelerationFunction( x_grid, y_grid, z_grid, ...
%     n_charges, n_masses);
% 
% 
% [xv2, accel] = velocityVerlet3D(ts, xv0(:,2), accelFunc(E_x, E_y, E_z));
% 
% 
% G = sqrt(hit_objective(xv).^2);
% 
% delta = 1e-6;
% 
% 
% for xx = 4:9
%     for yy = 4:9
%         fprintf('%i, %i\n', xx, yy);
% 
%         Ex2 = E_x;
%         Ex2(xx,yy,1) = Ex2(xx,yy,1) + delta;
%         [xv2, accel2] = velocityVerlet3D(ts, xv0, accelFunc(Ex2, E_y, E_z));
%         G2 = sqrt(hit_objective(xv2).^2);
% 
%         dGdEx_meas2(xx,yy) = dGdEx_meas2(xx,yy) + (G2-G)/delta;
% 
% 
% 
% 
%     end
% end
% 


% dGdEx_meas_sum = dGdEx_meas + dGdEx_meas2;



figure(40)
imagesc(x_grid, y_grid, dGdEx_meas');
axis xy image
colorbar
title('Meas') 
% figure(42)
% imagesc(x_grid, y_grid, dGdEx_meas2');
% axis xy image
% colorbar
% title('Meas') 
% figure(43)
% imagesc(x_grid, y_grid, (dGdEx_meas' + dGdEx_meas2'));
% axis xy image
% colorbar
% title('Meas') 
 %%

nParticle = 1;


hit_objective = @(x_v) hitObjective3D_wrap(...
            x_v, x_p, y_p, z_p, vx_p, vy_p, vz_p, obj_weights);
% hit_objective = @(x_v) hitObjective3D_wrap(...
%             x_v, x_p, y_p, z_p, vx_p, vy_p, vz_p, obj_weights);
V_Comsol = zeros(Nx, Ny, Nz, 3);
V_Comsol(:,:,:,1) = E_x; 
V_Comsol(:,:,:,2) = E_y; 
V_Comsol(:,:,:,3) = E_z;
%[dGdEx0, dGdEy0, dGdEz0, G0, xv00, DG, xv_dual] = VV_get_dual_E_final(n_charges, n_masses, E_x, E_y, E_z,  x_grid, y_grid, z_grid, xv0, nParticle, hit_objective);
[dGdEx_sum, dGdEy_sum, dGdEz_sum, G_sum, xv_all, DG_old, xv_dual, Nt, sysMatrix, iix2, iiy2, iiz2, w0002, w0012, w0102, w0112, w1002, w1012, w1102, w1112, S_p,  accelInterpMatrix_old, xv, DG_o, accelMatrix_old, dGdEx_old] ...
    = VV_get_dual_E_v21_objdebug(n_charges, n_masses, V_Comsol, ...
    x_grid, y_grid, z_grid, xv0, 1, hit_objective, ts);


figure(41)
imagesc(x_grid, y_grid, dGdEx_sum(:,:,1)');
axis xy image
colorbar
title('dual')
figure(50)
imagesc((dGdEx_sum(:,:,1)'-dGdEx_meas')./dGdEx_meas')
(dGdEx_sum(:,:,1)'-dGdEx_meas')./dGdEx_sum(:,:,1)'
%%