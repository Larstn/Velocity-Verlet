function [dF_dV, G, x_v_fw, dF_dE_x_d, dF_dE_y_d, dF_dE_z_d, dF_dE_x_fw] = solvewithComsolE_accuracy(x1_p, x2_p, x3_p,...
    nParticle, D_T, x1_init, x2_init, x3_init, v1_init, v2_init, v3_init,...
    E_x, E_y, E_z, m, q, x_grid, y_grid, z_grid)


[x_grid, y_grid, z_grid, d_x, d_y, d_z, N]     = Setup_Grid3DCom(x_grid, y_grid, z_grid);
[x_v_init]                                  = Setup_Particle3D(...
    x1_init,x2_init, x3_init, v1_init, v2_init, v3_init, nParticle);
%[E_x, E_y, E_z]                             = Setup_Fields3DGrad(x_grid,...
    %y_grid, z_grid, d_x, d_y, d_z, V);
[id_x1, id_x2, id_x3, id_v1, id_v2, id_v3]  = get_Index3D(nParticle);

% Interpolation
%
interp_E_x = @(x, y, z) interpn(x_grid, y_grid, z_grid, E_x, x, y, z,...
    'linear', 0);
interp_E_y = @(x, y, z) interpn(x_grid, y_grid, z_grid, E_y, x, y, z,...
    'linear', 0);
interp_E_z = @(x, y, z) interpn(x_grid, y_grid, z_grid, E_z, x, y, z,...
    'linear', 0);

%% Solve Systems
%
% Forward System
%
    [t_fw,x_v_fw] = forwardSystem_3D(interp_E_x, interp_E_y, interp_E_z,...
        x_v_init, D_T, nParticle, q, m);

% Solve Objective Function
%
    [G, D_G] = objectiveFunction_3D(x_v_fw, nParticle, x1_p, x2_p, x3_p);

% Solve Dual System 
%
    [t_d, xd] = dualSystem_3D(D_G, D_T, nParticle);
%%    
    dF_dE_x_fw = zeros(length(x_grid),length(y_grid),length(z_grid));


%% Calculate the measured Sensitivity 
%

%G_pert = zeros(1,6*nParticle);
G_t = sum(G(id_x1));
for ii = 1:length(x_grid)
    ii
    for jj = 1:length(y_grid)
        
        for kk = 1:length(z_grid)
            
            eps = 1e-3;
            
            
            E_x_pert = E_x;
            
            E_x_pert(ii,jj,kk) =  E_x(ii,jj,kk) + eps;
            
            interp_E_x_pert = @(x, y, z) interpn(x_grid, y_grid, z_grid,...
                E_x_pert, x, y, z, 'linear', 0);
            
            [t_Fw2, x_v_fw_pert] = forwardSystem_3D(interp_E_x_pert, ...
                interp_E_y, interp_E_z, x_v_init, D_T, nParticle, q, m);
            
            [G_pert,~] = objectiveFunction_3D(x_v_fw_pert,...
                nParticle, x1_p, x2_p, x3_p);
            G_pert_s = sum(G_pert(id_x1));
            dF_dE_x_fw(ii,jj,kk) = (G_pert_s - G_t) / eps;
            (G_pert_s - G_t) / eps;
            
            
            
            
        end
        
    end
end

% for ii = 1:N
%     
%     for jj = 1:N
%         
%         for kk = 1:N
%             
%             eps = 1e-3;
%             
%             
%             E_x_pert = E_x;
%             
%             E_x_pert(ii,jj,kk) =  E_x(ii,jj,kk) + eps;
%             
%             interp_E_x_pert = @(x, y, z) interpn(x_grid, y_grid, z_grid,...
%                 E_x_pert, x, y, z, 'linear', 0);
%             
%             [t_Fw2, x_v_fw_pert] = forwardSystem_3D(interp_E_x_pert, ...
%                 interp_E_y, interp_E_z, x_v_init, D_T, nParticle, q, m);
%             
%             [G_pert,~] = objectiveFunction_3D(x_v_fw_pert,...
%                 nParticle, x1_p, x2_p, x3_p);
%             G_pert_s = sum(G_pert(id_x1));
%             dF_dE_x_fw(ii,jj,kk) = (G_pert_s - G_t) / eps;
%             
%             
%             
%             
%         end
%         
%     end
% end

%% Interpolation in Time and Weight Recovery
%
x_d = spline(t_fw, x_v_fw(:,id_x1)', t_d);
y_d = spline(t_fw, x_v_fw(:,id_x2)', t_d);
z_d = spline(t_fw, x_v_fw(:,id_x3)', t_d);

[i_x, i_y, i_z, w000, w100, w010, w110, w001, w101, w011, w111] = ...
    trilinear_weights(x_d,y_d,z_d, x_grid, y_grid, z_grid, d_x, d_y, d_z, nParticle);

%% Get the Value of Integrant
%
d_t_s       = -centeredDiff(t_d);
dF_dE_x_t   = bsxfun(@times,d_t_s, xd(:,id_v1)*q/m);
dF_dE_y_t   = bsxfun(@times,d_t_s, xd(:,id_v2)*q/m);
dF_dE_z_t   = bsxfun(@times,d_t_s, xd(:,id_v3)*q/m);

%% Calculate Sensitivity with Dual Function 
%
[dF_dE_x_d]  = multipleRestriction3D(i_x, i_y, i_z,...
    w000, w100, w010,w110, w001, w101, w011, w111, ...
    dF_dE_x_t, size(E_x,1), size(E_x,2), size(E_x,3), nParticle );
[dF_dE_y_d]  = multipleRestriction3D(i_x, i_y, i_z,...
    w000, w100, w010,w110, w001, w101, w011, w111, ...
    dF_dE_y_t, size(E_x,1), size(E_x,2), size(E_x,3), nParticle );
[dF_dE_z_d]  = multipleRestriction3D(i_x, i_y, i_z,...
    w000, w100, w010,w110, w001, w101, w011, w111, ...
    dF_dE_z_t, size(E_x,1), size(E_x,2), size(E_x,3), nParticle );

%% Plots
%{
for i = 1:N
    
    figure()
    subplot(2,1,1)
    imagesc(x_grid, y_grid, dF_dE_x_fw(:,:,i)',[-0.1 0])
    colorbar
    
    hold on
    plot(x_v_fw(:, id_x1), x_v_fw(:, id_x2), 'wo-')
    title('Forward')
    
    subplot(2,1,2)
    imagesc(x_grid, y_grid, dF_dE_x_d(:,:,i)', [-0.1 0])
    colorbar
    
    hold on
    plot(x_v_fw(:, id_x1), x_v_fw(:, id_x2), 'wo-')
    title('Dual')
    
end
%}
%% Get Voltage Sensitivy
dF_dV = zeros(N(1),N(2),N(3));

dF_dV(3:end,:,:)   = dF_dV(3:end,:,:)     +...
    (-0.5/d_x)*dF_dE_x_d(2:end-1,:,:);
dF_dV(1:end-2,:,:) = dF_dV(1:end-2,:,:)   +...
    (0.5/d_x)*dF_dE_x_d(2:end-1,:,:);

dF_dV(:,3:end,:)   = dF_dV(:,3:end,:)     +...
    (-0.5/d_y)*dF_dE_y_d(:,2:end-1,:);
dF_dV(:,1:end-2,:) = dF_dV(:,1:end-2,:)   +...
    (0.5/d_y)*dF_dE_y_d(:,2:end-1,:);

dF_dV(:,:,3:end)   = dF_dV(:,:,3:end)     +...
    (-0.5/d_z)*dF_dE_z_d(:,:,2:end-1);
dF_dV(:,:,1:end-2) = dF_dV(:,:,1:end-2)   +...
    (0.5/d_z)*dF_dE_z_d(:,:,2:end-1);

end