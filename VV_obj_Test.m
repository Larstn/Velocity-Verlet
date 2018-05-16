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


ts = linspace(0, 1, Nt);



V           = repmat(peaks(Nx), [1 1 2]); assert(Nx == Ny);

V(:,:,2) = V(:,:,1);

E_x         = -centeredDiff(V, 1);
E_y         = -centeredDiff(V, 2);
E_z         = -centeredDiff(V, 3);

hit_objective = @(x_v) hitObjective3D_wrap(...
            x_v, x_p, y_p, z_p, vx_p, vy_p, vz_p, obj_weights);
        
        
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


figure(40)
imagesc(x_grid, y_grid, dGdEx_meas');
axis xy image
colorbar
title('Meas') 




nParticle = 1;


hit_objective = @(x_v) hitObjective3D_wrap(...
            x_v, x_p, y_p, z_p, vx_p, vy_p, vz_p, obj_weights);
        
        
        
        
        
        
        %%
particles = Particle(xv0(4:6,1),'Velocity',n_charges,n_masses,xv0(1:3,1), 0, 1, Nt);



%particles.accelerationFunction = accelFuncdebug(VV.Ex, VV.Ey, VV.Ez);

VV = Velocity_Verlet([x_grid(1), x_grid(end);y_grid(1),y_grid(end);z_grid(1),z_grid(end)], [Nx, Ny, Nz], 1, particles, hit_objective, E_x, E_y, E_z);

%accelFuncdebug = VV.accelerationFunction();
%accelFuncdebug = accelerationFunction_debug( x_grid, y_grid, z_grid); 

%VV.ParticleArray(1).accelerationFunction = accelFuncdebug(VV.Ex, VV.Ey, VV.Ez);


VV = VV.calculateF_dF;
disp('compare xv')
sum(VV.ParticleArray.xv == xv)
VV.plotdFdEx

figure(50)
%imagesc((VV.dFdEx(:,:,1)'-dGdEx_meas')./dGdEx_meas')
imagesc((VV.dFdEx(:,:,1)')./dGdEx_meas')

(VV.dFdEx(:,:,1)'-dGdEx_meas')./VV.dFdEx(:,:,1)'
