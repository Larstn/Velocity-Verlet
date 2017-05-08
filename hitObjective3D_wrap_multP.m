%% Calculation of G and dGdxv
% This function calculates the objective Function and its derivative with
% respect to the position/velocity vector (gradient)
function [G, dGdxv] = hitObjective3D_wrap_multP(xv, x_p, y_p, z_p,...
    vx_p, vy_p, vz_p, obj_weights)

    Nt = size(xv,1)/6;
    nParticle = size(xv,2);
    
    G_sum = 0;
    dGdxv_sum = zeros(6*Nt,1);
    
    for i = 1:nParticle
    
        [ix_x_relevant, ix_y_relevant, ix_z_relevant, iv_x_relevant,...
            iv_y_relevant, iv_z_relevant, relevant_pos] =...
            getrelevantidx(xv(:,i), x_p, Nt);

        [wobj_1, wobj_2] = objectiveInterpW(relevant_pos, x_p);

        [x_pos, y_pos, z_pos, x_vel, y_vel, z_vel] = getActualPos(xv,...
            wobj_1, wobj_2, ix_x_relevant, ix_y_relevant, ix_z_relevant, ...
            iv_x_relevant, iv_y_relevant, iv_z_relevant);

        [G, Gx, Gy, Gz, Gvx, Gvy, Gvz] = tHitObjective(x_pos, y_pos, z_pos,...
            x_vel, y_vel, z_vel, x_p, y_p, z_p, vx_p, vy_p, vz_p, obj_weights);

        [dw1dx1, dw1dx2, dw2dx1, dw2dx2] = weightDiff(relevant_pos, x_p);

        dGdxv = calculateDG(xv(:,i), x_p, y_p, z_p, vx_p, vy_p, vz_p, wobj_1, ...
            wobj_2, dw1dx1, dw2dx1, dw1dx2, dw2dx2, ix_x_relevant, ...
            ix_y_relevant, ix_z_relevant, iv_x_relevant, iv_y_relevant, ...
            iv_z_relevant, obj_weights);
    
        G = G_sum + G.^2;
        dGdxv_sum = dGdxv_sum + 1/nParticle*2*dGdxv*G;
    end
    %dGdxv2 = 2*G*dGdxv;


    G = sqrt(1/nParticle * G);
    
    
end