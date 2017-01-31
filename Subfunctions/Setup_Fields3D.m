function [E_x, E_y, E_z] = Setup_Fields3D(x_grid, y_grid, z_grid, d_x, d_y, d_z)



V1 = -repmat(y_grid, length(x_grid), 1);
%V = repmat(V, length(z_grid));

V = ones(length(x_grid), length(y_grid), length(z_grid));
for i = 1:length(z_grid)
    
    V(:,:,i) = V1;
    
end
%V = -ones(length(x_grid), length(y_grid), length(z_grid));

% Calculate the Electric Field

E_x = - centeredDiff(V, 1) / d_x;

E_y = - centeredDiff(V, 2) / d_y;

E_z = - centeredDiff(V, 3) / d_z;

end