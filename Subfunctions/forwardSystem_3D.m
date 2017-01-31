function [t_fw, x_v_fw] = forwardSystem_3D(interp_E_x, interp_E_y, interp_E_z, x_v_init, D_T, nParticle, q, m)

% Set up of the forward Problem

[id_x1, id_x2, id_x3, id_v1, id_v2, id_v3]= get_Index3D(nParticle);
%D_T_n = D_T(1):0.25:D_T(end);

D_T_n = [D_T(end) D_T(1)];
A = [zeros(3*nParticle), -eye(3*nParticle)
    zeros(3*nParticle), zeros(3*nParticle)];

f_fw = @(x_v) [zeros(nParticle,1);
               zeros(nParticle,1); 
               zeros(nParticle,1);
               (q/m)*interp_E_x(x_v(id_x1), x_v(id_x2), x_v(id_x3));
               (q/m)*interp_E_y(x_v(id_x1), x_v(id_x2), x_v(id_x3));
               (q/m)*interp_E_z(x_v(id_x1), x_v(id_x2), x_v(id_x3))];

function [x] = forward_Diff_eq(t, x_v_fw) 
    x = -A*x_v_fw + f_fw(x_v_fw);
   % keyboard
end

% Solving the Forward Problem with ODE45

[t_fw, x_v_fw] = ode45(@forward_Diff_eq, D_T_n, x_v_init);

end


