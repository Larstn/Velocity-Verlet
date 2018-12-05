% B_p
function [Bp] = B_matrix_primal_Array(ts, Dx_I, Dy_I, Dz_I, N_particles)

%Da_vec = [1 2 3 4 5 6 7 8 9 10];

%ts = 0:1:2;

Nt = length(ts);
dt = ts(2) - ts(1);

[ix_x, ix_y, ix_z, iv_x, iv_y, iv_z] = get_Index3D(Nt);
ia_x = ix_x;
ia_y = ix_y;
ia_z = ix_z;


% 
% Bp1_Ex = sparse([ix_x(3:Nt), ix_y(3:Nt), ix_z(3:Nt)],...
% [ix_x(2:Nt-1), ix_x(2:Nt-1), ix_x(2:Nt-1)],...
% [grad_Ex(ia_x(2:Nt-1)), grad_Ex(ia_y(2:Nt-1)), grad_Ex(ia_z(2:Nt-1))]*0.5*dt*dt, 6*Nt, 6*Nt);
% 
% Bp1_Ey = sparse([ix_x(3:Nt), ix_y(3:Nt), ix_z(3:Nt)],...
% [ix_y(2:Nt-1), ix_y(2:Nt-1), ix_y(2:Nt-1)],...
% [grad_Ey(ia_x(2:Nt-1)), grad_Ey(ia_y(2:Nt-1)), grad_Ey(ia_z(2:Nt-1))]*0.5*dt*dt, 6*Nt, 6*Nt);
% 
% Bp1_Ez = sparse([ix_x(3:Nt), ix_y(3:Nt), ix_z(3:Nt)],...
% [ix_z(2:Nt-1), ix_z(2:Nt-1), ix_z(2:Nt-1)],...
% [grad_Ez(ia_x(2:Nt-1)), grad_Ez(ia_y(2:Nt-1)), grad_Ez(ia_z(2:Nt-1))]*0.5*dt*dt, 6*Nt, 6*Nt);
% 
% 
% 
% Bp2_Ex = sparse([iv_x(3:Nt), iv_y(3:Nt), iv_z(3:Nt)],...
%     [ix_x(2:Nt-1), ix_x(2:Nt-1), ix_x(2:Nt-1)],...
%     [grad_Ex(ia_x(2:Nt-1)), grad_Ex(ia_y(2:Nt-1)), grad_Ex(ia_z(2:Nt-1))]*0.5*dt, 6*Nt, 6*Nt);
% 
% Bp2_Ey = sparse([iv_x(3:Nt), iv_y(3:Nt), iv_z(3:Nt)],...
%     [ix_y(2:Nt-1), ix_y(2:Nt-1), ix_y(2:Nt-1)],...
%     [grad_Ey(ia_x(2:Nt-1)), grad_Ey(ia_y(2:Nt-1)), grad_Ey(ia_z(2:Nt-1))]*0.5*dt, 6*Nt, 6*Nt);
% 
% Bp2_Ez = sparse([iv_x(3:Nt), iv_y(3:Nt), iv_z(3:Nt)],...
%     [ix_z(2:Nt-1), ix_z(2:Nt-1), ix_z(2:Nt-1)],...
%     [grad_Ex(ia_x(2:Nt-1)), grad_Ex(ia_y(2:Nt-1)), grad_Ex(ia_z(2:Nt-1))]*0.5*dt, 6*Nt, 6*Nt);
% 
% 
% 
% Bp3_Ex = sparse([iv_x(3:Nt), iv_y(3:Nt), iv_z(3:Nt)],...
%     [ix_x(3:Nt), ix_x(3:Nt), ix_x(3:Nt)],...
%     [grad_Ex(ia_x(3:Nt)), grad_Ex(ia_y(3:Nt)), grad_Ex(ia_z(3:Nt))]*0.5*dt, 6*Nt, 6*Nt);
% 
% Bp3_Ey = sparse([iv_x(3:Nt), iv_y(3:Nt), iv_z(3:Nt)],...
%     [ix_y(3:Nt), ix_y(3:Nt), ix_y(3:Nt)],...
%     [grad_Ey(ia_x(3:Nt)), grad_Ey(ia_y(3:Nt)), grad_Ey(ia_z(3:Nt))]*0.5*dt, 6*Nt, 6*Nt);
% 
% Bp3_Ez = sparse([iv_x(3:Nt), iv_y(3:Nt), iv_z(3:Nt)],...
%     [ix_z(3:Nt), ix_z(3:Nt), ix_z(3:Nt)],...
%     [grad_Ez(ia_x(3:Nt)), grad_Ez(ia_y(3:Nt)), grad_Ez(ia_z(3:Nt))]*0.5*dt, 6*Nt, 6*Nt);

idx1_vec = [];
idx2_vec = [];
val_vec = [];

for ii = 1:N_particles
    
    
    idx1 = [ix_x(3:Nt), ix_y(3:Nt), ix_z(3:Nt)] + (ii-1)*6*Nt;
    idx2 = [ix_x(2:Nt-1), ix_x(2:Nt-1), ix_x(2:Nt-1)] + (ii-1)*6*Nt;
    val = [Dx_I(ia_x(2:Nt-1)), Dx_I(ia_y(2:Nt-1)), Dx_I(ia_z(2:Nt-1))]*0.5*dt*dt;
    
    idx1_vec = [idx1_vec, idx1(:)'];
    idx2_vec = [idx2_vec, idx2(:)'];
    val_vec  = [val_vec, val(:)'];
    
    
% 
% Bp1_Ix = sparse([ix_x(3:Nt), ix_y(3:Nt), ix_z(3:Nt)],...
% [ix_x(2:Nt-1), ix_x(2:Nt-1), ix_x(2:Nt-1)],...
% [Dx_I(ia_x(2:Nt-1)), Dx_I(ia_y(2:Nt-1)), Dx_I(ia_z(2:Nt-1))]*0.5*dt*dt, 6*Nt, 6*Nt);

    idx1 = [ix_x(3:Nt), ix_y(3:Nt), ix_z(3:Nt)] + (ii-1)*6*Nt;
    idx2 = [ix_y(2:Nt-1), ix_y(2:Nt-1), ix_y(2:Nt-1)] + (ii-1)*6*Nt;
    val  = [Dy_I(ia_x(2:Nt-1)), Dy_I(ia_y(2:Nt-1)), Dy_I(ia_z(2:Nt-1))]*0.5*dt*dt;
    
    
    idx1_vec = [idx1_vec, idx1(:)'];
    idx2_vec = [idx2_vec, idx2(:)'];
    val_vec  = [val_vec, val(:)'];

% Bp1_Iy = sparse([ix_x(3:Nt), ix_y(3:Nt), ix_z(3:Nt)],...
% [ix_y(2:Nt-1), ix_y(2:Nt-1), ix_y(2:Nt-1)],...
% [Dy_I(ia_x(2:Nt-1)), Dy_I(ia_y(2:Nt-1)), Dy_I(ia_z(2:Nt-1))]*0.5*dt*dt, 6*Nt, 6*Nt);
% 

    idx1 = [ix_x(3:Nt), ix_y(3:Nt), ix_z(3:Nt)] + (ii-1)*6*Nt;
    idx2 = [ix_z(2:Nt-1), ix_z(2:Nt-1), ix_z(2:Nt-1)] + (ii-1)*6*Nt;
    val  = [Dz_I(ia_x(2:Nt-1)), Dz_I(ia_y(2:Nt-1)), Dz_I(ia_z(2:Nt-1))]*0.5*dt*dt;
    
    idx1_vec = [idx1_vec, idx1(:)'];
    idx2_vec = [idx2_vec, idx2(:)'];
    val_vec  = [val_vec, val(:)'];
    
% 
% Bp1_Iz = sparse([ix_x(3:Nt), ix_y(3:Nt), ix_z(3:Nt)],...
% [ix_z(2:Nt-1), ix_z(2:Nt-1), ix_z(2:Nt-1)],...
% [Dz_I(ia_x(2:Nt-1)), Dz_I(ia_y(2:Nt-1)), Dz_I(ia_z(2:Nt-1))]*0.5*dt*dt, 6*Nt, 6*Nt);
    

    idx1 = [iv_x(3:Nt), iv_y(3:Nt), iv_z(3:Nt)] + (ii-1)*6*Nt;
    idx2 = [ix_x(2:Nt-1), ix_x(2:Nt-1), ix_x(2:Nt-1)] + (ii-1)*6*Nt;
    val  = [Dx_I(ia_x(2:Nt-1)), Dx_I(ia_y(2:Nt-1)), Dx_I(ia_z(2:Nt-1))]*0.5*dt;
    
    idx1_vec = [idx1_vec, idx1(:)'];
    idx2_vec = [idx2_vec, idx2(:)'];
    val_vec  = [val_vec, val(:)'];    
% 
% Bp2_Ix = sparse([iv_x(3:Nt), iv_y(3:Nt), iv_z(3:Nt)],...
%     [ix_x(2:Nt-1), ix_x(2:Nt-1), ix_x(2:Nt-1)],...
%     [Dx_I(ia_x(2:Nt-1)), Dx_I(ia_y(2:Nt-1)), Dx_I(ia_z(2:Nt-1))]*0.5*dt, 6*Nt, 6*Nt);


    idx1 = [iv_x(3:Nt), iv_y(3:Nt), iv_z(3:Nt)] + (ii-1)*6*Nt;
    idx2 = [ix_y(2:Nt-1), ix_y(2:Nt-1), ix_y(2:Nt-1)] + (ii-1)*6*Nt;
    val  = [Dy_I(ia_x(2:Nt-1)), Dy_I(ia_y(2:Nt-1)), Dy_I(ia_z(2:Nt-1))]*0.5*dt;
    
    idx1_vec = [idx1_vec, idx1(:)'];
    idx2_vec = [idx2_vec, idx2(:)'];
    val_vec  = [val_vec, val(:)'];    
    
%     
% Bp2_Iy = sparse([iv_x(3:Nt), iv_y(3:Nt), iv_z(3:Nt)],...
%     [ix_y(2:Nt-1), ix_y(2:Nt-1), ix_y(2:Nt-1)],...
%     [Dy_I(ia_x(2:Nt-1)), Dy_I(ia_y(2:Nt-1)), Dy_I(ia_z(2:Nt-1))]*0.5*dt, 6*Nt, 6*Nt);


    idx1 = [iv_x(3:Nt), iv_y(3:Nt), iv_z(3:Nt)] + (ii-1)*6*Nt;
    idx2 = [ix_z(2:Nt-1), ix_z(2:Nt-1), ix_z(2:Nt-1)] + (ii-1)*6*Nt;
    val  = [Dz_I(ia_x(2:Nt-1)), Dz_I(ia_y(2:Nt-1)), Dz_I(ia_z(2:Nt-1))]*0.5*dt;
    
    idx1_vec = [idx1_vec, idx1(:)'];
    idx2_vec = [idx2_vec, idx2(:)'];
    val_vec  = [val_vec, val(:)']; 
    
% 
% Bp2_Iz = sparse([iv_x(3:Nt), iv_y(3:Nt), iv_z(3:Nt)],...
%     [ix_z(2:Nt-1), ix_z(2:Nt-1), ix_z(2:Nt-1)],...
%     [Dz_I(ia_x(2:Nt-1)), Dz_I(ia_y(2:Nt-1)), Dz_I(ia_z(2:Nt-1))]*0.5*dt, 6*Nt, 6*Nt);

    idx1 = [iv_x(3:Nt), iv_y(3:Nt), iv_z(3:Nt)] + (ii-1)*6*Nt;
    idx2 = [ix_x(3:Nt), ix_x(3:Nt), ix_x(3:Nt)] + (ii-1)*6*Nt;
    val  = [Dx_I(ia_x(3:Nt)), Dx_I(ia_y(3:Nt)), Dx_I(ia_z(3:Nt))]*0.5*dt;
    
    idx1_vec = [idx1_vec, idx1(:)'];
    idx2_vec = [idx2_vec, idx2(:)'];
    val_vec  = [val_vec, val(:)']; 

% 
% Bp3_Ix = sparse([iv_x(3:Nt), iv_y(3:Nt), iv_z(3:Nt)],...
%     [ix_x(3:Nt), ix_x(3:Nt), ix_x(3:Nt)],...
%     [Dx_I(ia_x(3:Nt)), Dx_I(ia_y(3:Nt)), Dx_I(ia_z(3:Nt))]*0.5*dt, 6*Nt, 6*Nt);

    idx1 = [iv_x(3:Nt), iv_y(3:Nt), iv_z(3:Nt)] + (ii-1)*6*Nt;
    idx2 = [ix_y(3:Nt), ix_y(3:Nt), ix_y(3:Nt)] + (ii-1)*6*Nt;
    val  = [Dy_I(ia_x(3:Nt)), Dy_I(ia_y(3:Nt)), Dy_I(ia_z(3:Nt))]*0.5*dt;
    
    idx1_vec = [idx1_vec, idx1(:)'];
    idx2_vec = [idx2_vec, idx2(:)'];
    val_vec  = [val_vec, val(:)']; 
 
% 
% Bp3_Iy = sparse([iv_x(3:Nt), iv_y(3:Nt), iv_z(3:Nt)],...
%     [ix_y(3:Nt), ix_y(3:Nt), ix_y(3:Nt)],...
%     [Dy_I(ia_x(3:Nt)), Dy_I(ia_y(3:Nt)), Dy_I(ia_z(3:Nt))]*0.5*dt, 6*Nt, 6*Nt);

    idx1 = [iv_x(3:Nt), iv_y(3:Nt), iv_z(3:Nt)] + (ii-1)*6*Nt;
    idx2 = [ix_z(3:Nt), ix_z(3:Nt), ix_z(3:Nt)] + (ii-1)*6*Nt;
    val  = [Dz_I(ia_x(3:Nt)), Dz_I(ia_y(3:Nt)), Dz_I(ia_z(3:Nt))]*0.5*dt;
    
    idx1_vec = [idx1_vec, idx1(:)'];
    idx2_vec = [idx2_vec, idx2(:)'];
    val_vec  = [val_vec, val(:)']; 

% 
% Bp3_Iz = sparse([iv_x(3:Nt), iv_y(3:Nt), iv_z(3:Nt)],...
%     [ix_z(3:Nt), ix_z(3:Nt), ix_z(3:Nt)],...
%     [Dz_I(ia_x(3:Nt)), Dz_I(ia_y(3:Nt)), Dz_I(ia_z(3:Nt))]*0.5*dt, 6*Nt, 6*Nt);



    ia_x = ia_x + 3*Nt;
    ia_y = ia_y + 3*Nt;
    ia_z = ia_z + 3*Nt;

end

Bp = sparse(idx1_vec, idx2_vec, val_vec,ii*6*Nt, ii*6*Nt);

%Bp = Bp1_Ex + Bp2_Ex + Bp3_Ex + Bp1_Ey + Bp2_Ey + Bp3_Ey + Bp1_Ez + Bp2_Ez + Bp3_Ez + Bp1_Ix + Bp2_Ix + Bp3_Ix + Bp1_Iy + Bp2_Iy + Bp3_Iy + Bp1_Iz + Bp2_Iz + Bp3_Iz;

end