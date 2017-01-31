%% Get relevant indicies that give the index where the target was hit
% Uses the target position and returns the corresponding indices that are
% the first occasion where the particle hits the plane descriped by the
% x-postion of the target. If the target gets not hit, it returns the index
% that corresponds to the closest position. 
%% 
function  [ix_x_relevant, ix_y_relevant, ix_z_relevant, iv_x_relevant, iv_y_relevant, iv_z_relevant, relevant_pos]= getrelevantidx(xv, x_p, Nt)

[ix_x, ix_y, ix_z, iv_x, iv_y, iv_z] = get_Index3D(Nt);

dist_x = abs(xv(ix_x) - x_p);
[dist_xsort, idx_sort] = sort(dist_x);
idxdiff = abs(idx_sort - idx_sort(1));
last_idx = find(idxdiff == 1);

%always choose the first time the trajectory hits the x wall
relevant_idx = sort(idx_sort(1:last_idx));
relevant_pos = xv(relevant_idx(1:2));

ix_x_relevant = ix_x(relevant_idx(1:2));
ix_y_relevant = ix_x_relevant + Nt;
ix_z_relevant = ix_y_relevant + Nt;
iv_x_relevant = ix_z_relevant + Nt;
iv_y_relevant = iv_x_relevant + Nt;
iv_z_relevant = iv_y_relevant + Nt;
end