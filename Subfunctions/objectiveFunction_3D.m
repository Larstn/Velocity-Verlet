function [G, D_G] = objectiveFunction_3D(x_v_fw, nParticle, x1_p, x2_p, x3_p)

[id_x1, id_x2, id_x3, id_v1, id_v2, id_v3] = get_Index3D(nParticle);

%x_p = 0.3;
%y_p = 1;

%vx_p = 0;
%vy_p = 0;

G = zeros(1,6*nParticle);
D_G = zeros(1,6*nParticle);

G(id_x1)=  (x1_p - x_v_fw(end,id_x1)).^2; % (y_p - x_v_fw(end,2)).^2;
G(id_x2)=  (x2_p - x_v_fw(end,id_x2)).^2; % (y_p - x_v_fw(end,2)).^2;
G(id_x3)=  (x3_p - x_v_fw(end,id_x3)).^2; % (y_p - x_v_fw(end,2)).^2;

D_G(id_x1) = 2*(x1_p - x_v_fw(end,id_x1)); % 2*(y_p - x_v_fw(end,2))];
D_G(id_x2) = 2*(x2_p - x_v_fw(end,id_x2));
D_G(id_x3) = 2*(x3_p - x_v_fw(end,id_x3));


end