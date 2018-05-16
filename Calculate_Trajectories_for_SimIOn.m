E_x = zeros(581,1401,2);
E_y = zeros(581,1401,2);
E_z = zeros(581,1401,2);




% E_x(:,:,1) = Ex;
% E_x(:,:,2) = Ex;
% 
% E_y(:,:,1) = Ey;
% E_y(:,:,2) = Ey;
% 
% E_z(:,:,1) = Ez; 
% E_z(:,:,2) = Ez;

E_x(:,:,1) = Exfull5_old;
E_x(:,:,2) = Exfull5_old;

E_y(:,:,1) = Eyfull5_old;
E_y(:,:,2) = Eyfull5_old;

E_z(:,:,1) = Ezfull5_old; 
E_z(:,:,2) = Ezfull5_old;

n_charges = 1;
n_masses = 100;
KE = 125;
electron_mass       = 1.6605e-27;
elementary_charge   = 1.60217662e-19;

speed = sqrt(2*KE*elementary_charge/(n_masses*electron_mass));%VVVx(1)*1e3;%sqrt(2*KE*elementary_charge/electron_mass);
xinit = 0*ones(1,11);
yinit = linspace(-5e-3,5e-3 ,11);
zinit = 0*ones(1,11);
vxinit = speed*ones(1,11);
vyinit = 0*ones(1,11);
vzinit = 0*ones(1,11);
xv0 = [xinit; yinit; zinit; vxinit; vyinit; vzinit];
nParticle = 11;

% x_grid = x3(1:2901);
% x_grid = x_grid';
% y_grid = y3(1:2901:end);
% y_grid = y_grid';
x_grid = x16(1:581);
x_grid = x_grid';
y_grid = y16(1:581:end);
y_grid = y_grid';
z_grid = [-0.001 0.001];

ts = linspace(0,19.2e-6,300);
obj_weights = [1,1,1,0,0,0];
x1_p = x3(2901)*ones(1,11);
x2_p = 0*x1_p;
x3_p = 0*x1_p;
v1_p = 0*x1_p;
v2_p = 0*x1_p;
v3_p = 0*x1_p;
hit_objective = @(xv) hitObjective3D_wrap(xv, x1_p, x2_p, x3_p, v1_p, v2_p, v3_p, obj_weights);


[dGdEx_sum, dGdEy_sum, dGdEz_sum, G_sum, xv_all, DG, xv_dual, Nt] ...
    = VV_get_dual_E_v20(n_charges, n_masses, E_x, E_y, E_z, ...
    x_grid, y_grid, z_grid, xv0, nParticle, hit_objective, ts);
figure(101)
imagesc(E_x(:,:,1))
figure(102)
imagesc(E_y(:,:,1))


[ix_x, ix_y, ix_z, id_v1, id_v2, id_v3, ~] =...
                     get_Index3D(Nt);

figure(100)
hold all
for i = 1:nParticle
    plot(ts,xv_all(id_v1,i))
end
% figure(2)
% plot(xv_all(id_v2))
% figure(3)
% plot(xv_all(id_v3))
                 
                 
figure(310)
clf
hold on 
for i = 1:11
    plot3(xv_all(ix_x,i), xv_all(ix_y,i), xv_all(ix_z,i))
end

xlabel('x [mm]')
ylabel('y [mm]')
title('Trajectories from Comsol/VV')

figure(1002)
clf(1002)
hold on
imagesc(y_grid, x_grid,E_x(:,:,1))
for i = 1:11
    plot(xv_all(ix_y(1:end),i), xv_all(ix_x(1:end),i),'ro')
    
end

xlabel('x [mm]')
ylabel('y [mm]')
title('Trajectories from Comsol/VV')



