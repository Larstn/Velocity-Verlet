%% Calculation of the Derivative of the G vector
% Calculates the derivative of G towards the input vector with respect to 
% all the components of the input vector (gradient of the input vector)

function DGdxv = calculateDG(xv, x_p, y_p, z_p,...
    vx_p, vy_p, vz_p, wobj_1, wobj_2, dw1dx1, dw2dx1, dw1dx2, dw2dx2, ...
    ix_x_relevant, ix_y_relevant, ix_z_relevant, iv_x_relevant,...
    iv_y_relevant, iv_z_relevant, obj_weights)

w_y = obj_weights(2);
w_z = obj_weights(3);
w_vx = obj_weights(4);
w_vy = obj_weights(5);
w_vz = obj_weights(6);

DGdxv = xv*0;

DGdxv(ix_x_relevant(1)) = DGdx(xv, x_p, y_p, z_p,...
    vx_p, vy_p, vz_p, wobj_1, wobj_2, dw1dx1, dw2dx1, ...
    ix_x_relevant, ix_y_relevant, ix_z_relevant, iv_x_relevant,...
    iv_y_relevant, iv_z_relevant, obj_weights);
DGdxv(ix_x_relevant(2)) = DGdx(xv, x_p, y_p, z_p,...
    vx_p, vy_p, vz_p, wobj_1, wobj_2, dw1dx2, dw2dx2, ...
    ix_x_relevant, ix_y_relevant, ix_z_relevant, iv_x_relevant,...
    iv_y_relevant, iv_z_relevant, obj_weights);

[ DGdxv(ix_y_relevant(1)), DGdxv(ix_y_relevant(2))]  = ...
    DGdxv_23(w_y ,ix_y_relevant, wobj_1, wobj_2, xv, y_p);
[ DGdxv(ix_z_relevant(1)), DGdxv(ix_z_relevant(2))]  = ...
    DGdxv_23(w_z ,ix_z_relevant, wobj_1, wobj_2, xv, z_p);
[ DGdxv(iv_x_relevant(1)), DGdxv(iv_x_relevant(2))]  = ...
    DGdxv_23(w_vx ,iv_x_relevant, wobj_1, wobj_2, xv, vx_p);
[ DGdxv(iv_y_relevant(1)), DGdxv(iv_y_relevant(2))]  = ...
    DGdxv_23(w_vy ,iv_y_relevant, wobj_1, wobj_2, xv, vy_p);
[ DGdxv(iv_x_relevant(1)), DGdxv(iv_x_relevant(2))]  = ...
    DGdxv_23(w_vz ,iv_z_relevant, wobj_1, wobj_2, xv, vz_p);

end