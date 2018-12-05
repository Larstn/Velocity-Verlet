function [G_sum, dGdxv] = hitObjective3D_pP_wrap(xv, x_p, y_p, z_p,...
    vx_p, vy_p, vz_p, obj_weights, N_particles)


    G_sum = 0;
    dGdxv = zeros(size(xv,1), N_particles);

    for ii = 1:N_particles
        
        [G, dGdxv(:,ii)] = hitObjective3D_wrap(xv(:,ii), x_p, y_p, z_p,...
    vx_p, vy_p, vz_p, obj_weights);

        G_sum = G_sum + G^2;
        dGdxv(:,ii) = dGdxv(:,ii) * 2 * 1/N_particles * G;
        
        
    end 

        G_sum = sqrt(1/N_particles * G_sum);
        
        dGdxv = dGdxv .* 0.5 .* (1/G_sum);



end 