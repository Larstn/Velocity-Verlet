% how would i like to use the VV code

% initial conditions

Nparticles = 1234567;

xyz = [ ... ] % column vectors?
%or
xyz = zeros(Nparticles, 3);
velocity_xyz = [ ... ]  % column vectors?
energies = [ ... ]

% Might want to run like....

vv.simulate(xyz, velocity_xyz, potentials)
%or
vv.simulate(xyz, velocity_xyz, electricField)

% or
vv.simulate(..., dt = 1e-6)  or (relative_dt = 1e-3) % find some smart criterion


% Extract trajectories?

trajectories = vv.getAnswer() % ?

or 

trajectories = vv.simulate( ... )

% need to think about why one or the other is better
% do i want results to be stored in the VV object/module or just
% sent to the user and forgotten?