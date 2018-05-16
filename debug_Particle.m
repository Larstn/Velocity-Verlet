particles = Particle(xv0(4:6,1),'Velocity',n_charges,n_masses,xv0(1:3,1), 0, 1, Nt);


accelFuncdebug = accelerationFunction_debug( x_grid, y_grid, z_grid); 


 particles.accelerationFunction = accelFuncdebug(VV.Ex, VV.Ey, VV.Ez);
 particles = particles.runTrajectory();

 %particles.xv - xv
 
 sum(particles.xv == xv)
 
 
 particles2 = VV.ParticleArray(1);
 partilces2.accelerationFunction = accelFuncdebug(VV.Ex, VV.Ey, VV.Ez);
 partilces2 = particles.runTrajectory();
 
 sum(partilces2.xv == xv)
 
 %%
 
 [iix, iiy, iiz, w000, w001, w010, w011, w100, w101, w110, w111] ...
                    = VV.trilinear_weights(1);
                
sum(iix == iix2)
sum(iiy == iiy2)
sum(iiz == iiz2)

sum(w000 == w0002)
sum(w001 == w0012)
sum(w010 == w0102)
sum(w011 == w0112)
sum(w100 == w1002)
sum(w101 == w1012)
sum(w110 == w1102)
sum(w111 == w1112)


max(max(S_p - VV.S_p))
%%
xv_dual - VV.ParticleArray(1).xv_dual
%%

isequal(xv_all,VV.ParticleArray(1).xv)

%%
accelInterpMatrix2 = VV.get_accelInterpmatrix(iix, iiy, iiz, ...
                    w000, w001, w010, w011, w100, w101, w110, w111, ...
                    1);
%VV.x_grid - x_grid

accelInterpMatrix2 - accelInterpMatrix

%%
VV.DG - DG

%%

[G2, DG2] = VV.Fobj(xv_all);

sum(VV.DG - DG2)

sum(VV.DG - hit_objective(xv_all))

%%
clc
VV.Fobj = hit_objective;
[G, DG] = VV.Fobj(xv_all);
[G2, DG2] = hit_objective(xv_all);
sum(DG - DG2)
sum(DG - DG_old)
%%
clc
VV = Velocity_Verlet([x_grid(1), x_grid(end);y_grid(1),y_grid(end);z_grid(1),z_grid(end)], [Nx, Ny, Nz], 1, particles, hit_objective, E_x, E_y, E_z);
[G, DG] = VV.Fobj(xv_all);
[G2, DG2] = hit_objective(xv_all);
sum(DG - DG2)

VV = VV.calculateF_dF;
sum(DG2 - VV.DG)

%%
sum(VV.DG - DG_o)
%%
sum(DG_old - hit_objective(xv_all))
%%

sum(VV.accelMatrix(:) - accelMatrix_old(:))
sum(VV.accelInterpMatrix(:) - accelInterpMatrix_old(:))

%%

sum(dGdEx_old(:) - VV.dGdEx(:))


%% 


sum(VV.dGdEx_sum(:) - dGdEx_sum(:))


