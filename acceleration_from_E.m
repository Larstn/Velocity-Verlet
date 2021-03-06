function accelFunc = accelerationFunction( x_grid, y_grid, z_grid)

interpolant = @(E, x, y, z) interpn(x_grid, y_grid, z_grid, E, x, y, z, 'linear', 0);


accelFunc = @(Ex,Ey,Ez) @(t, xyz) [...
    interpolant(Ex, xyz(1), xyz(2), xyz(3)); ...
    interpolant(Ey, xyz(1), xyz(2), xyz(3)); ...
    interpolant(Ez, xyz(1), xyz(2), xyz(3)); ...
    ]*(elementary_charge/electron_mass);

end