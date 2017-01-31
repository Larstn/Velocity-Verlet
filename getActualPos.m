%% Return the actual position 
% Returns the actual position of the particle when it hits the target plane
%% 

function [x_pos, y_pos, z_pos, x_vel, y_vel, z_vel] = ...
    getActualPos(xv, wobj_1, wobj_2, ...
    ix_x_relevant, ix_y_relevant, ix_z_relevant,...
    iv_x_relevant, iv_y_relevant, iv_z_relevant)

    x_pos = wobj_1*xv(ix_x_relevant(1)) + wobj_2*xv(ix_x_relevant(2));
    y_pos = wobj_1*xv(ix_y_relevant(1)) + wobj_2*xv(ix_y_relevant(2));
    z_pos = wobj_1*xv(ix_z_relevant(1)) + wobj_2*xv(ix_z_relevant(2));

    x_vel = wobj_1*xv(iv_x_relevant(1)) + wobj_2*xv(iv_x_relevant(2));
    y_vel = wobj_1*xv(iv_y_relevant(1)) + wobj_2*xv(iv_y_relevant(2));
    z_vel = wobj_1*xv(iv_z_relevant(1)) + wobj_2*xv(iv_z_relevant(2));
    
end