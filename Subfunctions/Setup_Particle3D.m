function [x_v_init] = Setup_Particle3D(x1_init, x2_init, x3_init, v1_init, v2_init, v3_init, nParticle)

% Initial Values for Speed and Position

[id_x1, id_x2, id_x3, id_v1, id_v2, id_v3] = get_Index3D(nParticle);

x_v_init = zeros(6*nParticle, 1);

x_v_init(id_x1) = x1_init; %x-position
x_v_init(id_x2) = x2_init; %y-position
x_v_init(id_x3) = x3_init; %z-position

x_v_init(id_v1) = v1_init; %v_x
x_v_init(id_v2) = v2_init; %v_y
x_v_init(id_v3) = v3_init; %v_z

end
