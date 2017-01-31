function [dGdEx_sum, dGdEy_sum, dGdEz_sum, G_sum, xv_all, DG, xv_dual] = VV_get_dual_E_mult_units(n_charges, n_masses, E_x, E_y, E_z, x_grid, y_grid, z_grid, xv0, x_p, y_p, z_p, vx_p, vy_p, vz_p, nParticle, obj_weights)
    elementary_charge   = 1.60217662e-19;
    electron_mass       = 9.1093856e-31;

    Nx = size(E_x, 1); % No units
    Ny = size(E_x, 2); % No units
    Nz = size(E_y, 3);
    interpolant = @(E, x, y, z) interpn(x_grid, y_grid, z_grid, E, x, y, z, 'linear', 0);
    accelFunc = @(Ex,Ey,Ez) @(t, xyz) [...
        interpolant(Ex, xyz(1), xyz(2), xyz(3)); ...
        interpolant(Ey, xyz(1), xyz(2), xyz(3)); ...
        interpolant(Ez, xyz(1), xyz(2), xyz(3)); ...
        ]*((n_charges*elementary_charge)/(n_masses*electron_mass)); % -electroncharge for negative charge [m*s^-2], thus electron charge / electron mass

    Nt = 100;
    ts = linspace(0, 1, Nt); % second

    [ix_x, ix_y, ix_z, iv_x, iv_y, iv_z, ia_x, ia_y, ia_z] = get_Index3D(Nt);


    G_sum = 0;
    dGdEx_sum = 0*E_x;
    dGdEy_sum = 0*E_y;
    dGdEz_sum = 0*E_z;
    xv_all    = zeros(6*Nt, nParticle);

    for ii = 1:nParticle
        
        xv_start = xv0(:,ii);
        %disp('in get Dual')
        %disp(xv_start)
        %[systemMatrix, initMatrix, accelMatrix] = velocityVerletMatrices3D(ts);
        
        [xv, accel] = velocityVerlet3D(...
            ts, xv_start, accelFunc(E_x, E_y, E_z));

        %xv = recreate_units(xv, 1, 1, 0);

        [iix, iiy, iiz, w000, w001, w010, w011, w100, w101, w110, w111] ...
            = trilinear_weights(xv(ix_x), xv(ix_y), xv(ix_z), ...
                x_grid, y_grid, z_grid);

        sz = [Nx Ny Nz];
        i000 = sub2ind(sz, iix, iiy, iiz);
        i001 = sub2ind(sz, iix+1, iiy, iiz);
        i010 = sub2ind(sz, iix, iiy+1, iiz);
        i011 = sub2ind(sz, iix+1, iiy+1, iiz);
        i100 = sub2ind(sz, iix, iiy, iiz+1);
        i101 = sub2ind(sz, iix+1, iiy, iiz+1);
        i110 = sub2ind(sz, iix, iiy+1, iiz+1);    
        i111 = sub2ind(sz, iix+1, iiy+1, iiz+1);

        nn = 1:Nt;
        accelInterpMatrix = sparse([nn, nn, nn, nn, nn, nn, nn, nn], ...
            [i000; i001; i010; i011; i100; i101; i110; i111], ...
            [w000; w001; w010; w011; w100; w101; w110; w111], ...
            Nt, numel(E_x));
% 
        accelInterpMatrix = get_accelInterpmatrix(iix, iiy, iiz, ...
            w000, w001, w010, w011, w100, w101, w110, w111, ...
            Nx, Ny, Nz, Nt);

        [Ix, Iy, Iz] = get_I(xv, ix_x, ix_y, ix_z,...
            x_grid, y_grid, z_grid, Nx, Ny, Nz, Nt);
        %ax = accelInterpMatrix * E_x(:);
        %ay = accelInterpMatrix * E_y(:);
        %az = accelInterpMatrix * E_z(:);

%         Ix = trilinear_weights_I_x(xv(ix_x), xv(ix_y), xv(ix_z), ...
%             x_grid, y_grid, z_grid, Nx, Ny, Nz, Nt);
%         Iy = trilinear_weights_I_y(xv(ix_x), xv(ix_y), xv(ix_z), ...
%             x_grid, y_grid, z_grid, Nx, Ny, Nz, Nt);
%         Iz = trilinear_weights_I_z(xv(ix_x), xv(ix_y), xv(ix_z), ...
%             x_grid, y_grid, z_grid, Nx, Ny, Nz, Nt);
%         %xv_new = systemMatrix \ (-initMatrix*xv_start - accelMatrix*[ax; ay; az]);

        [G, DG] = hitObjective3D_wrap(xv, x_p, y_p, z_p, vx_p, vy_p, vz_p, obj_weights);
        
        [systemMatrix, initMatrix, accelMatrix] = velocityVerletMatrices3D(ts);
        %systemMatrix = systemMatrix*(1/((n_charges*elementary_charge)/(n_masses*electron_mass)));
         systemMatrix = systemMatrix *(1/((n_charges*elementary_charge)/(n_masses*electron_mass)));

        
        [S_p] = get_PrimalS(Ix, Iy, Iz, E_x, E_y, E_z, Nt, ts, systemMatrix);
        
%         D_I_1          = [Ix*E_x(:); Ix*E_y(:); Ix*E_z(:)]*(elementary_charge/electron_mass); % [kg*m*s^-2*A^-1]
%         D_I_2          = [Iy*E_x(:); Iy*E_y(:); Iy*E_z(:)]*(elementary_charge/electron_mass);
%         D_I_3          = [Iz*E_x(:); Iz*E_y(:); Iz*E_z(:)]*(elementary_charge/electron_mass);
% 
%         grad_E_1 = zeros(1,3*Nt);
%         grad_E_2 = zeros(1,3*Nt);
%         grad_E_3 = zeros(1,3*Nt);
%     
%         Bp = B_matrix_primal(ts,grad_E_1, grad_E_2, grad_E_3, D_I_1, D_I_2, D_I_3);
%         S_p = systemMatrix + Bp;
        %S_p = S_p * (1/(elementary_charge/electron_mass));

        
        xv_dual = S_p' \ DG;

        dGda = xv_dual' * (-accelMatrix);
        dGdax = dGda(ix_x);
        dGdEx = reshape(dGdax * accelInterpMatrix, [Nx Ny Nz]);

        dGday = dGda(ix_y);
        dGdEy = reshape(dGday * accelInterpMatrix, [Nx Ny Nz]);

        dGdaz = dGda(ix_z);
        dGdEz = reshape(dGdaz * accelInterpMatrix, [Nx Ny Nz]);

        G_sum = G_sum + G;
        dGdEx_sum = dGdEx_sum + dGdEx;
        dGdEy_sum = dGdEy_sum + dGdEy;
        dGdEz_sum = dGdEz_sum + dGdEz;
        xv_all(:,ii) = xv;
        %xv_all(1,ii)

    end

end

