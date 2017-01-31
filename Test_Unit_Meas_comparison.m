%% Input parameters for this test case

Nx          = 11;
Ny          = 11;
Nz          = 2;

Nt = 100;

xv0 = [0; 0; 0; .8; .499; 0]*(elementary_charge/electron_mass);

x_p = 0.25*(elementary_charge/electron_mass);
y_p  = 0.3*(elementary_charge/electron_mass);
z_p  = 0*(elementary_charge/electron_mass);
vx_p = 1*(elementary_charge/electron_mass);
vy_p = 0*(elementary_charge/electron_mass);
vz_p = 0*(elementary_charge/electron_mass);

obj_weights = [1, 1, 1, 0, 0, 0];

%% 

%%
% Setting natural constants


elementary_charge   = 1.60217662e-19;
electron_mass       = 9.1093856e-31;
%%


ts = linspace(0, 1, Nt);


x_grid = linspace(-1, 1, Nx)*(elementary_charge/electron_mass);
y_grid = linspace(-1, 1, Ny)*(elementary_charge/electron_mass);
z_grid = linspace(-1, 1, Nz)*(elementary_charge/electron_mass);

%V           = ones(Nx, Ny, Nz);
%V           = 3*randn(Nx, Ny, Nz);
V           = repmat(peaks(Nx), [1 1 2]); assert(Nx == Ny);
V(:,:,2) = V(:,:,1);
E_x         = -centeredDiff(V, 1);
E_y         = -centeredDiff(V, 2);
E_z         = -centeredDiff(V, 3);

% put these out here so I can use them...
interpolant = @(E, x, y, z) interpn(x_grid, y_grid, z_grid, E, x, y, z, 'linear', 0);
accelFunc = @(Ex,Ey,Ez) @(t, xyz) [...
    interpolant(Ex, xyz(1), xyz(2), xyz(3)); ...
    interpolant(Ey, xyz(1), xyz(2), xyz(3)); ...
    interpolant(Ez, xyz(1), xyz(2), xyz(3)); ...
    ]*(elementary_charge/electron_mass);


[ix_x, ix_y, ix_z, iv_x, iv_y, iv_z, ia_x, ia_y, ia_z] = get_Index3D(Nt);



[xv, accel] = velocityVerlet3D(ts, xv0, accelFunc(E_x, E_y, E_z));



[G, dGdxv] = hitObjective3D(xv, obj_weights);


[systemMatrix, initMatrix, accelMatrix] = velocityVerletMatrices3D(ts);


ax1 = accel(:,1);
ay1 = accel(:,2);
az1 = accel(:,3);

delta = 1e-9;
dGdEx_meas = 0*E_x(:,:,1);

%for xx = 1:Nx
%for yy = 1:Ny
for xx = 5:8
for yy = 5:8
    fprintf('%i, %i\n', xx, yy);
    
    Ex2 = E_x;
    Ex2(xx,yy,1) = Ex2(xx,yy,1) + delta;
    [xv2, accel2] = velocityVerlet3D(ts, xv0, accelFunc(Ex2, E_y, E_z));
    [G2] = hitObjective3D(xv2, obj_weights);
%     Ex2 = E_x;
%     Ex2(xx,yy,1) = Ex2(xx,yy,1) - delta;
%     [xv2] = velocityVerlet3D(ts, xv0, accelFunc(Ex2, E_y, E_z));
%     [G3] = sumObjective3D(xv2);
    
    dGdEx_meas(xx,yy) = (G2-G)/delta;
    
%     ax2 = accelInterpMatrix * Ex2(:);
%     ay2 = accelInterpMatrix * E_y(:);
%     az2 = accelInterpMatrix * E_z(:);
    
%     figure(1); clf
%     plot(ax1);
%     hold on
%     plot(ax2);
%     legend('Orig', 'Perturbed')
    %pause

%     figure(10); clf
%     imagesc(x_grid, y_grid, Ex2(:,:,1)');
%     axis xy image vis3d

%     colorbar
%     ax = axis;
%     hold on
%     plot(xv2(ix_x), xv2(ix_y), 'b-');
%     pause(0.01)
    
    
end
end

figure(4)
subplot(3,1,2)
imagesc(x_grid, y_grid, dGdEx_meas');
axis xy image
colorbar
title('Meas')


nParticle = 1;


[dGdEx0, dGdEy0, dGdEz0, G0, xv00, DG, xv_dual] = VV_get_dual_E_mult_units(E_x, E_y, E_z,  x_grid, y_grid, z_grid, xv0, x_p, y_p, z_p, vx_p, vy_p, vz_p, nParticle, obj_weights);


figure(4)
subplot(3,1,1)
imagesc(x_grid, y_grid, dGdEx0(:,:,1)');
axis xy image
colorbar
title('dual')


Nx          = 11;
Ny          = 11;
Nz          = 2;



V           = repmat(peaks(Nx), [1 1 2]); assert(Nx == Ny);
V(:,:,2) = V(:,:,1);
E_x         = -centeredDiff(V, 1);
E_y         = -centeredDiff(V, 2);
E_z         = -centeredDiff(V, 3);
T_span      = [0 1];

%nParticle = 1;
%xv0 = zeros(6,nParticle);
%xv0 = [0; 0; 0; .8; .499; 0];



[dGdEx1, dGdEy, dGdEz, G, xv1, xv_new] = VV_get_dual_E_mult(E_x, E_y, E_z,  x_grid, y_grid, z_grid, xv0, x_p, y_p, z_p, vx_p, vy_p, vz_p, nParticle, obj_weights);

figure(4)
subplot(3,1,3)
imagesc(x_grid, y_grid, dGdEx1(:,:,1)');
axis xy image
colorbar
title('dual')