%% Creation of the acceleration Function
% Creates the acceleration function with the given grid as well as the
% charge (as a multiple of the elementary charge) and the masss (as
% multiple of the electron mass). The created acceleration function uses
% the electric field in all three components as well as the position as an
% input
%%
function accelFunc = accelerationFunction_debug( x_grid, y_grid, z_grid)


%elementary_charge   = 1.60217662e-19;
%electron_mass       = 1.6605e-27;


interpolant = @(E, x, y, z) interpn(x_grid, y_grid, z_grid, E, x, y, z, 'linear', 0);


accelFunc = @(Ex,Ey,Ez) @(t, xyz) [...
    interpolant(Ex, xyz(1), xyz(2), xyz(3)); ...
    interpolant(Ey, xyz(1), xyz(2), xyz(3)); ...
    interpolant(Ez, xyz(1), xyz(2), xyz(3)); ...
    ];

end