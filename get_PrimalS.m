%% Calculation of the Primal Matrix
% This function calculates the matrix of the Primal System (for a
% derivation with respect to the electric field).
function [S_p] = get_PrimalS(Ix, Iy, Iz, E_x, E_y, E_z, Nt, ts, systemMatrix)
%     elementary_charge   = 1.60217662e-19;
%     electron_mass       = 9.1093856e-31;
%     n_charges = 1;
%     n_masses = 1;
    D_I_1          = [Ix*E_x(:); Ix*E_y(:); Ix*E_z(:)]; % [kg*m*s^-2*A^-1]
    D_I_2          = [Iy*E_x(:); Iy*E_y(:); Iy*E_z(:)];
    D_I_3          = [Iz*E_x(:); Iz*E_y(:); Iz*E_z(:)];

    grad_E_1 = zeros(1,3*Nt);
    grad_E_2 = zeros(1,3*Nt);
    grad_E_3 = zeros(1,3*Nt);

    Bp = B_matrix_primal(ts,grad_E_1, grad_E_2, grad_E_3, D_I_1, D_I_2, D_I_3);
    

    S_p = systemMatrix + Bp;
end