%% Calculation of G and dGdxv
% This function calculates the objective Function and its derivative with
% respect to the position/velocity vector (gradient)
function [G, dGdxv] = hitObjective3D_wrap(xv, x_p, y_p, z_p,...
    vx_p, vy_p, vz_p, obj_weights)

    Nt = size(xv,1)/6;

    [ix_x_relevant, ix_y_relevant, ix_z_relevant, iv_x_relevant,...
        iv_y_relevant, iv_z_relevant, relevant_pos] =...
        getrelevantidx(xv, x_p, Nt);
    if xv(ix_x_relevant(2)) < x_p
        wobj_1 = 0;
        wobj_2 = 1;
        
    else
    [wobj_1, wobj_2] = objectiveInterpW(relevant_pos, x_p);
    
    end

    [x_pos, y_pos, z_pos, x_vel, y_vel, z_vel] = getActualPos(xv,...
        wobj_1, wobj_2, ix_x_relevant, ix_y_relevant, ix_z_relevant, ...
        iv_x_relevant, iv_y_relevant, iv_z_relevant);

    [G, Gx, Gy, Gz, Gvx, Gvy, Gvz] = tHitObjective(x_pos, y_pos, z_pos,...
        x_vel, y_vel, z_vel, x_p, y_p, z_p, vx_p, vy_p, vz_p, obj_weights);
    
    if xv(ix_x_relevant(2)) < x_p
        w_x  = obj_weights(1);
        w_y  = obj_weights(2);
        w_z  = obj_weights(3);
        w_vx = obj_weights(4);
        w_vy = obj_weights(5);
        w_vz = obj_weights(6);
        dGdxv = 0 * xv;
        dGdxv(ix_x_relevant(2)) = 2*w_x*(x_pos - x_p);
        dGdxv(ix_y_relevant(2)) = 2*w_y*(y_pos - y_p);
        dGdxv(ix_z_relevant(2)) = 2*w_z*(z_pos - z_p);
        dGdxv(iv_x_relevant(2)) = 2*w_vx*(x_vel - vx_p);
        dGdxv(iv_y_relevant(2)) = 2*w_vy*(y_vel - vy_p);
        dGdxv(iv_z_relevant(2)) = 2*w_vz*(z_vel - vz_p);
    else 
        [dw1dx1, dw1dx2, dw2dx1, dw2dx2] = weightDiff(relevant_pos, x_p);
    

        dGdxv = calculateDG(xv, x_p, y_p, z_p, vx_p, vy_p, vz_p, wobj_1, ...
            wobj_2, dw1dx1, dw2dx1, dw1dx2, dw2dx2, ix_x_relevant, ...
            ix_y_relevant, ix_z_relevant, iv_x_relevant, iv_y_relevant, ...
            iv_z_relevant, obj_weights);

    end
end