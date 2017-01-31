%% Hit Objective 3D with Target Input
% Returns the value of the objective function G as well as its
% differentiation towards all components of the xv-vector
%%
function [G, dGdxv] = hitObjective3Dtarget(xv, obj_weights, x_p, y_p,...
    z_p, vx_p, vy_p, vz_p)

Nt = size(xv,1)/6;

[ix_x_relevant, ix_y_relevant, ix_z_relevant, iv_x_relevant,...
    iv_y_relevant, iv_z_relevant, relevant_pos]= ...
    getrelevantidx(xv, x_p, Nt);


[wobj_1, wobj_2] = objectiveInterpW(relevant_pos, x_p);

[x_pos, y_pos, z_pos, x_vel, y_vel, z_vel] = getActualPos(xv,...
    wobj_1, wobj_2, ix_x_relevant, ix_y_relevant, ix_z_relevant, ...
    iv_x_relevant, iv_y_relevant, iv_z_relevant);

[G, Gx, Gy, Gz, Gvx, Gvy, Gvz] = tHitObjective(x_pos, y_pos, z_pos,...
    x_vel, y_vel, z_vel, x_p, y_p, z_p, vx_p, vy_p, vz_p, obj_weights);


[dw1dx1, dw1dx2, dw2dx1, dw2dx2] = weightDiff(relevant_pos, x_p);

% 
% [dx_posdx1, dy_posdx1, dz_posdx1, dx_veldx1, dy_veldx1, dz_veldx1] = getActualPos(xv,...
%     dw1dx1, dw2dx1, ix_x_relevant, ix_y_relevant, ix_z_relevant, ...
%     iv_x_relevant, iv_y_relevant, iv_z_relevant);
% 
% 
% [dx_posdx2, dy_posdx2, dz_posdx2, dx_veldx2, dy_veldx2, dz_veldx2] = getActualPos(xv,...
%     dw1dx2, dw2dx2, ix_x_relevant, ix_y_relevant, ix_z_relevant, ...
%     iv_x_relevant, iv_y_relevant, iv_z_relevant);

dGdxv = calculateDG(xv, x_p, y_p, z_p,...
    vx_p, vy_p, vz_p, wobj_1, wobj_2, dw1dx1, dw2dx1, dw1dx2, dw2dx2,...
    ix_x_relevant, ix_y_relevant, ix_z_relevant, iv_x_relevant,...
    iv_y_relevant, iv_z_relevant, obj_weights);


end