%% Input parameters for this test case
clc 
clear all
%%
% Setting natural constants
% in case that they are needed

elementary_charge   = 1.60217662e-19;
%electron_mass       = 9.1093856e-31;
atomic_mass = 1.6605e-27;

Nx          = 11;
Ny          = 11;
Nz          = 2;

Nt = 100;

x_grid = linspace(-1, 1, Nx)*(elementary_charge/atomic_mass);
y_grid = linspace(0, 1, round(0.5*Ny))*(elementary_charge/atomic_mass);
z_grid = linspace(-1, 1, Nz)*(elementary_charge/atomic_mass);

xv0 = [0.2 0.2; 0.2 0.2; -1 -1;0 0;0 0; 0 0]*(elementary_charge/atomic_mass);

x_p = 0.8*(elementary_charge/atomic_mass);
y_p  = -0.2*(elementary_charge/atomic_mass);
z_p  = 0;
vx_p = 1;
vy_p = 0;
vz_p = 0;

obj_weights = [1, 1, 0, 0, 0, 0];

n_charges = 1;
n_masses = 1;
nParticle = 2;

d_x = x_grid(2) - x_grid(1);
d_y = y_grid(2) - y_grid(1);
d_z = z_grid(2) - z_grid(1);

%% 
% Creating test environment with given Input

ts = linspace(0, 5e3, Nt);

% 
% V = repmat(linspace(5,-5,11)',1,6,2);
% V = V + repmat(linspace(-2.5,2.5,6),11,1,2);
% 
% V = [V(:,end:-1:2,:) V];
% %V           = repmat(peaks(Nx), [1 1 2]); assert(Nx == Ny);
% r_grid = [-y_grid(end:-1:2) y_grid];
% %V(:,:,2) = V(:,:,1);
% 
% E_x         = -centeredDiff(V, 1);
% E_r         = -centeredDiff(V, 2);
% E_z         = -centeredDiff(V, 3);

E_x = ones(Nx,Ny,Nz);
E_r = -ones(Nx,Ny,Nz)*0.7;
E_z = zeros(Nx,Ny,Nz);
V = -cumsum(E_x,1) - cumsum(E_r,2) - cumsum(E_z,3);


E_x         = -centeredDiff(V, 1) / d_x;
E_r         = -centeredDiff(V, 2) / d_y;
E_z         = -centeredDiff(V, 3) / d_z;
% E_2r = [E_r(:,end:-1:1,:) E_r(:,:,:)];
% r2_grid = [r_grid(end:-1:1) r_grid];

%%
    figure(20)
    subplot(3,1,1)
    imagesc(E_x(:,:,1))
    subplot(3,1,2)
    imagesc(E_r(:,:,1))
    subplot(3,1,3)
    imagesc(E_z(:,:,1))
    
%%
    nParticle = 2;


hit_objective = @(x_v) hitObjective3D_wrap(...
            x_v, x_p, y_p, z_p, vx_p, vy_p, vz_p, obj_weights);
% hit_objective = @(x_v) hitObjective3D_wrap(...
%             x_v, x_p, y_p, z_p, vx_p, vy_p, vz_p, obj_weights);
%V_Comsol = zeros(Nx, Ny, Nz, 3);
V_Comsol = cell(3,1);
V_Comsol{1} = E_x;
V_Comsol{2} = E_r; 
V_Comsol{3} = E_z;
%[dGdEx0, dGdEy0, dGdEz0, G0, xv00, DG, xv_dual] = VV_get_dual_E_final(n_charges, n_masses, E_x, E_y, E_z,  x_grid, y_grid, z_grid, xv0, nParticle, hit_objective);
[dGdEx_sum, dGdEy_sum, dGdEz_sum, dGdV, dGdV_xr, G_sum, xv_all, DG, xv_dual, Nt] ...
    = VV_get_dual_E_cylindrical(n_charges, n_masses, V_Comsol, ...
    x_grid, r_grid, z_grid, xv0, nParticle, hit_objective, ts);
%%

figure(401)
clf
imagesc(x_grid, r_grid, dGdEx_sum(:,:,1)');
axis xy image
colorbar
title('dual Ex')
hold on
plot(x_p,y_p,'rx')
plot(xv_all(1:100,1),xv_all(101:200,1),'k')
plot(xv_all(1:100,2),xv_all(101:200,2),'b')
figure(402)
clf
imagesc(x_grid, r_grid, dGdEy_sum(:,:,1)');
axis xy image
colorbar
title('dual Ey')
hold on
plot(x_p,y_p,'rx')
plot(xv_all(1:100,1),xv_all(101:200,1),'k')
plot(xv_all(1:100,2),xv_all(101:200,2),'b')
%%    

hit_objective = @(x_v) hitObjective3D_wrap(...
            x_v, x_p, y_p, z_p, vx_p, vy_p, vz_p, obj_weights);


%%
nParticle = 2;
dGdEx_meas = 0*E_x(:,:,1);

accelFunc = accelerationFunction( x_grid, r_grid, z_grid, ...
    n_charges, n_masses);


[xv, accel] = velocityVerlet3D(ts, xv0(:,1), accelFunc(E_x, E_r, E_z));
[xv2, accel] = velocityVerlet3D(ts, xv0(:,2), accelFunc(E_x, E_r, E_z));


G = sqrt(1/nParticle*(hit_objective(xv).^2 + hit_objective(xv2).^2));
G
delta = 1e-3;


for xx = 3:10%1:length(x_grid)
    for yy = 3:10%1:length(y_grid)
        fprintf('%i, %i\n', xx, yy);

        Ex2 = E_x;
        Ex2(xx,yy,1) = Ex2(xx,yy,1) + delta*E_x(xx,yy,1);
        max(max(Ex2-E_x))
        [xv3, accel2] = velocityVerlet3D(ts, xv0, accelFunc(Ex2, E_r, E_z));
        [xv4, accel2] = velocityVerlet3D(ts, xv0, accelFunc(Ex2, E_r, E_z));

        G2 = sqrt(1/nParticle*(hit_objective(xv3).^2 + hit_objective(xv4).^2));
        G2
        G-G2
        max(xv-xv3)
        max(xv-xv4)
        dGdEx_meas(xx,yy) = dGdEx_meas(xx,yy) + (G2-G)/(delta*E_x(xx,yy,1));




    end
end
%%

figure(44)
imagesc(dGdEx_meas')
axis xy image

figure(40)
clf
%imagesc(dGdEx_meas')
imagesc(x_grid, r_grid, dGdEx_meas');
axis xy image
hold on
plot(x_p,y_p,'rx')
plot(xv(1:100,1),xv(101:200,1),'k')
plot(xv3(1:100,1),xv3(101:200,1),'b')
%axis xy image
colorbar
title('Meas') 


figure(500)
imagesc((dGdEx_sum(:,:,1)'-dGdEx_meas')./dGdEx_meas')
colorbar
%%
dGdEr_meas = 0*E_r(:,:,1);


for xx = 3:10%1:length(x_grid)
    for yy = 3:10%1:length(y_grid)
        fprintf('%i, %i\n', xx, yy);

        Er2 = E_r;
        Er2(xx,yy,1) = Er2(xx,yy,1) + delta*E_r(xx,yy,1);
        max(max(Er2-E_x))
        [xv3, accel2] = velocityVerlet3D(ts, xv0, accelFunc(E_x, Er2, E_z));
        [xv4, accel2] = velocityVerlet3D(ts, xv0, accelFunc(E_x, Er2, E_z));

        G2 = sqrt(1/nParticle*(hit_objective(xv3).^2 + hit_objective(xv4).^2));
        G2
        G-G2
        max(xv-xv3)
        max(xv-xv4)
        dGdEr_meas(xx,yy) = dGdEr_meas(xx,yy) + (G2-G)/(delta*E_r(xx,yy,1));




    end
end

%%
figure(607)
imagesc(dGdEr_meas')
axis xy image

figure(608)
imagesc(dGdEy_sum(:,:,1)')
axis xy image

figure(609)
imagesc(dGdEr_meas'./dGdEy_sum(:,:,1)')
axis xy image
colorbar
%%
dGdEz_meas = 0*E_z(:,:,1);


for xx = 3:10%1:length(x_grid)
    for yy = 3:10%1:length(y_grid)
        fprintf('%i, %i\n', xx, yy);

        Ez2 = E_z;
        Ez2(xx,yy,1) = Ez2(xx,yy,1) + delta*E_z(xx,yy,1);
        max(max(Ez2-E_z))
        [xv3, accel2] = velocityVerlet3D(ts, xv0, accelFunc(E_x, E_r, Ez2));
        [xv4, accel2] = velocityVerlet3D(ts, xv0, accelFunc(E_x, E_r, Ez2));

        G2 = sqrt(1/nParticle*(hit_objective(xv3).^2 + hit_objective(xv4).^2));
        G2
        G-G2
        max(xv-xv3)
        max(xv-xv4)
        dGdEz_meas(xx,yy) = dGdEz_meas(xx,yy) + (G2-G)/(delta*E_z(xx,yy,1));




    end
end

%%
figure(607)
imagesc(dGdEz_meas')
axis xy image

figure(608)
imagesc(dGdEz_sum(:,:,1)')
axis xy image

figure(609)
imagesc(dGdEz_meas'./dGdEz_sum(:,:,1)')
axis xy image
colorbar




%%
delta = 1e-3;
dGdV_meas = 0*V(:,:,1);


for xx = 1:11%1:length(x_grid)
    for yy = 1:11%1:length(y_grid)
          fprintf('%i, %i\n', xx, yy);
            UV2 = V;
            
            accelFunc = accelerationFunction( x_grid, r_grid, z_grid, ...
    n_charges, n_masses);

            UV2(xx,yy,1) = UV2(xx,yy,1) + delta;
            max(max(UV2(:,:,1)-V(:,:,1)))
            E_x2         = -centeredDiff(UV2, 1) / d_x;
            E_r2         = -centeredDiff(UV2, 2) / d_y;
            E_z2         = -centeredDiff(UV2, 3) / d_z;
            
            [xv3, accel2] = velocityVerlet3D(ts, xv0(:,1), accelFunc(E_x2, E_r2, E_z2));
            [xv4, accel2] = velocityVerlet3D(ts, xv0(:,2), accelFunc(E_x2, E_r2, E_z2));
            
            G2 = sqrt(1/nParticle*(hit_objective(xv3).^2 + hit_objective(xv4).^2));
            G2
            G-G2
            max(xv-xv3)
            max(xv2-xv4)
            dGdV_meas(xx,yy) = dGdV_meas(xx,yy) + (G2-G)/(delta);
            
            
    end
end

%%

figure(6070)
imagesc(dGdV_meas')
axis xy image

figure(608)
imagesc(dGdV(:,:,1)')
axis xy image

figure(609)
imagesc(dGdV_meas'./dGdV(:,:,1)')
axis xy image
colorbar

%%
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


%%

% dGdEx_sum./dGdEx_meas
% (dGdEx_sum(:,:,1)'-dGdEx_meas')./dGdEx_sum(:,:,1)'
figure(501)
imagesc(x_grid,r_grid,dGdV(:,:,1)')
xlabel('x')
ylabel('r')
figure(502)
imagesc(x_grid,y_grid,dGdV_xr(:,:,1)')
%%