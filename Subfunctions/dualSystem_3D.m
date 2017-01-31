function [t_d, xd] = dualSystem_3D(D_G, D_T, nParticle)

[id_x1, id_x2, id_x3, id_v1, id_v2, id_v3]= get_Index3D(nParticle);

%D_T_n = D_T(1):0.25:D_T(end);
D_T_n = [D_T(end) D_T(1)];
% Transpose Matrix
A_d = [zeros(3*nParticle), zeros(3*nParticle);
    -eye(3*nParticle), zeros(3*nParticle)];

% g depends on the Objective Function
g = [zeros(nParticle,1);
     zeros(nParticle,1);
     zeros(nParticle,1);
     zeros(nParticle,1);
     zeros(nParticle,1);
     zeros(nParticle,1)];

dual_Diff_eq = @(t, xv)  g - A_d*xv;

%D_T_d = D_T(end:-1:1);
[t_d_reversed, x_d_reversed] = ode45(dual_Diff_eq, D_T_n, D_G);

t_d = D_T(2)- (t_d_reversed(end:-1:1) - D_T(1));
%t_d = D_T(2)- (t_d_reversed(end:-1:1) - D_T(1));

xd = x_d_reversed(end:-1:1,:);