function [i_x, i_y, i_z, w000, w100, w010, w110, w001, w101, w011, w111] = trilinear_weights(x,y,z, x_grid, y_grid, z_grid)
% trilinear_weights    Indices and weights for trilinear interpolation
%
% trilinear_weights(x, y, z, xgrid, ygrid, zgrid)
% Calculate weights and indices needed to interpolate a function defined
% on a regular 3D grid to nParticle lists of coordinates.
%
% Inputs:
%    x, y, z            N-element
% 
% Outputs:
%   i_x, i_y, i_z       Indices to evaluate function at
%   w000, ..., w111     Weights to assign to function values

assert(isvector(x))
assert(isvector(y))
assert(isvector(z))
assert(isvector(x_grid))
assert(isvector(y_grid))
assert(isvector(z_grid))
assert(isequal(size(x), size(y), size(z)))

% Coerce x_grid etc. to be row vectors if x, y, z are row vectors...
if isrow(x) && ~isrow(x_grid)
    row = @(A) reshape(A, 1, []);
    x_grid = row(x_grid);
    y_grid = row(y_grid);
    z_grid = row(z_grid);
elseif iscolumn(x) && ~iscolumn(x_grid)
    col = @(A) reshape(A, [], 1);
    x_grid = col(x_grid);
    y_grid = col(y_grid);
    z_grid = col(z_grid);
end


d_x = x_grid(2)-x_grid(1);
d_y = y_grid(2)-y_grid(1);
d_z = z_grid(2)-z_grid(1);

% TODO: what does "log" mean in "log_vec"? - It was supposed to mean
% "Logic" just because it is a binary vector
log_vec = ((x >= x_grid(1)) & (x <= x_grid(end)) & (y >= y_grid(1)) & (y <= y_grid(end)) & (z >= z_grid(1)) & (z <= z_grid(end)));

i_x = ones(size(x));
i_y = ones(size(y));
i_z = ones(size(z));

i_x(log_vec) = floor( (x(log_vec) - x_grid(1))/d_x) +1;
i_y(log_vec) = floor( (y(log_vec) - y_grid(1))/d_y) +1;
i_z(log_vec) = floor( (z(log_vec) - z_grid(1))/d_z) +1;

% Handle a special case: when x == x_grid(end) we round its position DOWN
% instead of UP, to allow all boundary values to be defined as one might
% expect.  Now x == x_grid(1) is in the first cell AND x == x_grid(end)
% is in the last cell.

i_x(x == x_grid(end)) = length(x_grid)-1;
i_y(y == y_grid(end)) = length(y_grid)-1;
i_z(z == z_grid(end)) = length(z_grid)-1;

x_c = x_grid(i_x);
y_c = y_grid(i_y);
z_c = z_grid(i_z);

% Recover weights

wx = ( x - x_c ) ./ d_x;
wy = ( y - y_c ) ./ d_y;
wz = ( z - z_c ) ./ d_z; 

w000 = (1-wx).*(1-wy).*(1-wz);
w100 = wx.*(1-wy).*(1-wz);
w010 = (1-wx).*wy.*(1-wz);
w110 = wx.*wy.*(1-wz);
w001 = (1-wx).*(1-wy).*wz;
w101 = wx.*(1-wy).*wz;
w011 = (1-wx).*wy.*wz;
w111 = wx.*wy.*wz;

w000(~log_vec) = 0;
w100(~log_vec) = 0;
w010(~log_vec) = 0;
w110(~log_vec) = 0;

w001(~log_vec) = 0;
w101(~log_vec) = 0;
w011(~log_vec) = 0;
w111(~log_vec) = 0;

assert(isequal(size(i_x), size(x), size(i_y), size(y), size(i_z), size(z), ...
    size(w000), size(w001), size(w010), size(w011), ...
    size(w100), size(w101), size(w110), size(w111)));

end