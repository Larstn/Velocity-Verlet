%% Calculate Spot Size

function [z_spot_fin, y_spot_fin, rms_r_spot] = calculateSpotsize(x_v_fw, Nt, Ex, Ey, Ez, xs, ys, zs, k_max,ts, n_charges, n_masses, nParticle, x_p, y_p, z_p)


[id_x1, id_x2, id_x3, id_v1, id_v2, id_v3] = get_Index3D(Nt);

t_y = x_v_fw(id_x2,:);
t_z = x_v_fw(id_x3,:); 

y_var = var(t_y, 0,2);

[y_spot, idx_y ]= min(y_var);

z_var = var(t_z,0,2);

[z_spot, idx_z] = min(z_var);


r = sqrt((t_y-y_p).^2 + (t_z-z_p).^2);

r_var = var(r,0,2);

[r_spot, idx_r] = min(r);


min_r = min(idx_r);
max_r = max(idx_r); 

relevant_radii = r(min_r:max_r, :);

rms_relevant_radii = rms(relevant_radii,2);

[rms_r_spot, idx_rms_r] = min(rms_relevant_radii);

idx_rms_r = idx_rms_r + min(idx_r) - 1;

ts_n = ts;
xv_n = x_v_fw;
[ix_xr, ix_yr, ix_zr, iv_xr, iv_yr, iv_zr] =...
                     get_Index3D(Nt);

                 
                 
                 
for k = 1:k_max
    N_r = (max_r - min_r)*10;
    td = 0;
    ts_n = linspace(ts_n(min_r-td),ts_n(max_r+td),N_r);

    xv_start = [xv_n(ix_xr(min_r-td),:); xv_n(ix_yr(min_r-td),:); xv_n(ix_zr(min_r-td),:); xv_n(iv_xr(min_r-td),:); xv_n(iv_yr(min_r-td),:); xv_n(iv_zr(min_r-td),:)]; 
    
    
    accelFunc = accelerationFunction( xs, ys, zs, ...
        n_charges, n_masses);
    xv_n= zeros(6*N_r,nParticle);
    for ii = 1:nParticle
        [xv_n(:,ii), ~] = velocityVerlet3D(...
                    ts_n, xv_start(:,ii), accelFunc(Ex, Ey, Ez));

    end

    [ix_xr, ix_yr, ix_zr, iv_xr, iv_yr, iv_zr] =...
                         get_Index3D(N_r);
                     
                     
    figure(55)
    clf
    hold on 
    for ll = 1:nParticle
        plot3(xv_n(ix_xr,ll), xv_n(ix_yr, ll), xv_n(ix_zr,ll))
    end 
    view(3)
    xlabel('x [m]')
    ylabel('y [m]')
    zlabel('z [m]')
    title('Spot Zoom')
    
    [id_x1, id_x2, id_x3, id_v1, id_v2, id_v3] = get_Index3D(Nt);

    t_y = xv_n(ix_yr,:);
    t_z = xv_n(ix_zr,:); 

    y_var = var(t_y, 0,2);

    [y_spot, idx_y ]= min(y_var);

    z_var = var(t_z,0,2);

    [z_spot, idx_z] = min(z_var);


    r = sqrt((t_y-y_p).^2 + (t_z-z_p).^2);

    r_var = var(r,0,2);

    [r_spot, idx_r] = min(r);


    min_r = min(idx_r);
    max_r = max(idx_r); 

    relevant_radii = r(min_r:max_r, :);

    rms_relevant_radii = rms(relevant_radii,2);

    [rms_r_spot, idx_rms_r] = min(rms_relevant_radii);

    idx_rms_r = idx_rms_r + min(idx_r) - 1;
    title_str = sprintf('Spot Zoom, RMS = %f mm', rms_r_spot*1e3);
    title(title_str)

    y_var_fin = var(t_y, 0,2);

    [y_spot_fin, idx_y ]= min(y_var_fin);

    z_var_fin = var(t_z,0,2);

    [z_spot_fin, idx_z] = min(z_var_fin);

    
end 

%%
% ts_n = ts;
% xv_n = x_v_fw;
% [ix_xr, ix_yr, ix_zr, iv_xr, iv_yr, iv_zr] =...
%                      get_Index3D(Nt);
% 
% for k = 1:k_max
% 
% 
% N_r = 20;
% td = 1;
% ts_n = linspace(ts_n(idx_z-td),ts_n(idx_z+td),N_r);
% xv_start = [xv_n(ix_xr(idx_z-td),:); xv_n(ix_yr(idx_z-td),:); xv_n(ix_zr(idx_z-td),:); xv_n(iv_xr(idx_z-td),:); xv_n(iv_yr(idx_z-1),:); xv_n(iv_zr(idx_z-1),:)]; 
%            
%     
% accelFunc = accelerationFunction( xs, ys, zs, ...
%     n_charges, n_masses);
% xv_n= zeros(6*N_r,nParticle);
% for ii = 1:nParticle
% [xv_n(:,ii), ~] = velocityVerlet3D(...
%                 ts_n, xv_start(:,ii), accelFunc(Ex, Ey, Ez));
%             
% end
% 
% [ix_xr, ix_yr, ix_zr, iv_xr, iv_yr, iv_zr] =...
%                      get_Index3D(N_r);
%                  
% %%
%                  
% figure(55)
% clf
% hold on 
% for ll = 1:nParticle
%     plot3(xv_n(ix_xr,ll), xv_n(ix_yr, ll), xv_n(ix_zr,ll))
% end 
% view(3)
% xlabel('x [m]')
% ylabel('y [m]')
% zlabel('z [m]')
% title('Spot Zoom')
% %%
% 
% 
% t_y = xv_n(ix_yr,:);
% t_z = xv_n(ix_zr,:); 
% 
% y_var_fin = var(t_y, 0,2);
% 
% [y_spot_fin, idx_y ]= min(y_var_fin);
% 
% z_var_fin = var(t_z,0,2);
% 
% [z_spot_fin, idx_z] = min(z_var_fin);
% 
% %pause(0.5)
% end
% %%
% % [outDFx2, outDFy2, outDFz2, outF2, x_v_fw2, DG2, xv_dual2, Nt2] = VV_get_dual_E_v23(n_charges, n_masses, V_Comsol, xs, ys, zs, xv_start, nParticle, hit_objective, ts_n); 
% % 
% % figure(56)
% % clf
% hold on 
% for ll = 1:nParticle
%     plot3(xv_n(ix_xr,ll), xv_n(ix_yr, ll), xv_n(ix_zr,ll))
% end 
% view(3)
% xlabel('x [m]')
% ylabel('y [m]')
% zlabel('z [m]')