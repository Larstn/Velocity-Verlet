function [x_grid, y_grid, z_grid, d_x, d_y, d_z, N] = Setup_Grid3DCom(x_grid, y_grid, z_grid)
% Number of Data Points

% Set up the Gridpoints for Interpolation etc.
%x_grid = linspace(-1, 1, N);
%y_grid = linspace(-1, 1, N);
%z_grid = linspace(-1, 1, N);

% Calculate the size of each box
d_x = x_grid(end) - x_grid(end-1);
d_y = y_grid(end) - y_grid(end-1);
d_z = z_grid(end) - z_grid(end-1);

N = [length(x_grid), length(y_grid), length(z_grid)];

end