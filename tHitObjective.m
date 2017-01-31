%% Hit Objective Values
% Returns the values of the Hit objective function for each direction as
% well as the overall sum
%%
function [G, Gx, Gy, Gz, Gvx, Gvy, Gvz] = tHitObjective(x_pos, y_pos, z_pos, x_vel, y_vel, z_vel, x_p, y_p, z_p, vx_p, vy_p, vz_p, obj_weights)

w_x  = obj_weights(1);
w_y  = obj_weights(2);
w_z  = obj_weights(3);
w_vx = obj_weights(4);
w_vy = obj_weights(5);
w_vz = obj_weights(6);

Gx = w_x*(x_pos - x_p).^2;
Gy = w_y*(y_pos - y_p).^2;
Gz = w_z*(z_pos - z_p).^2;
Gvx = w_vx*(x_vel - vx_p).^2;
Gvy = w_vy*(y_vel - vy_p).^2;
Gvz = w_vz*(z_vel - vz_p).^2;

G = sum([Gx Gy Gz Gvx Gvy Gvz]);


end 