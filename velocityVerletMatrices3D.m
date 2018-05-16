%% Velocity Verlet Matrix Creation
% calculates all the velocity verlet matrices
%%
function [systemMatrix, A0, B] = velocityVerletMatrices3D(ts)
% velocityVerletMatrices3D    Matrices and indexing for 3D velocity Verlet method
%
% systemMatrix, A0, B = velocityVerletMatrices3D(ts)
%
% Let xv be the position and velocity vector we seek.  Then
%
% systemMatrix*xv + initialMatrix*[x0; y0; z0; vx0; vy0; vz0] + accelMatrix*accelerations
%
% can be solved for xv.  Here,
%   length(ts) = Nt + 1
%   size(xv) = [6*Nt, 1]
%   size(accelerations) = [3*(Nt+1), 1]     i.e. [ax0; ... ; axNt; ay0; ... ... azNt]
%
% Use the indexing functions to get out the positions and velocities, i.e.
%   x(n) = xv(i_x(n))
%   v(n) = xv(i_v(n))
%

%% Establish an indexing convention
Nt = length(ts);
dt = ts(2) - ts(1);

[ix_x, ix_y, ix_z, iv_x, iv_y, iv_z] = get_Index3D(Nt);

ia_x = ix_x;
ia_y = ix_y;
ia_z = ix_z;

%% System matrix for 1D (x) system
% We'll turn it into the 3D matrix eventually.
% Look at velocityVerletMatrices1D for explanation...

unos = @(n) ones(1,n);

% ------- Diagonal
% It's a big old identity matrix.

ii_diag = [ix_x(1:Nt), iv_x(1:Nt)];
jj_diag = [ix_x(1:Nt), iv_x(1:Nt)];
vv_diag = [-unos(Nt), -unos(Nt)];

% ------- A
% A = [1, dt; 0, 1]


ii_A = [ix_x(2:Nt), ix_x(2:Nt), iv_x(2:Nt)];
jj_A = [ix_x(1:Nt-1), iv_x(1:Nt-1), iv_x(1:Nt-1)];
vv_A = [unos(Nt-1), dt*unos(Nt-1), unos(Nt-1)];


% ------- Build the system matrix!

numRows = 6*Nt;
numCols = 6*Nt;

systemMatrix = sparse(...
    [ii_A, ii_diag, ii_A+Nt, ii_diag+Nt, ii_A+2*Nt, ii_diag+2*Nt], ...
    [jj_A, jj_diag, jj_A+Nt, jj_diag+Nt, jj_A+2*Nt, jj_diag+2*Nt], ...
    [vv_A, vv_diag, vv_A, vv_diag, vv_A, vv_diag],...
    numRows, numCols);
%spy(systemMatrix)


%% Right-hand side

% ------- B
% B = [0.5*dt^2, 0; 0.5*dt, 0.5*dt]
% The B matrix is 2*Nt x Nt+1.  Yeah it's weird-shaped.  It turns Nt+1
% accelerations into Nt timestep updates.

ii_B = [ix_x(2:Nt), iv_x(2:Nt), iv_x(2:Nt)];
jj_B = [ia_x(1:Nt-1), ia_x(1:Nt-1), ia_x(2:Nt)];
vv_B = [0.5*dt*dt*unos(Nt-1), 0.5*dt*unos(Nt-1), 0.5*dt*unos(Nt-1)];

numRows = 6*Nt;
numCols = 3*Nt;
B = sparse([ii_B, ii_B+Nt, ii_B+2*Nt],...
    [jj_B, jj_B+Nt, jj_B+2*Nt], ...
    [vv_B, vv_B, vv_B], numRows, numCols);
%figure(3); clf
%imagesc(B)

% ------- Initial conditions
% The first timestep needs to take initial conditions:
% [x(1); v(1)] = A*[x(0); v(0)] + B*[a(0); a(1)]
% So we need an RHS vector of the right size, 2*Nt x 1.
% Let's get RHS = A0 * [x(0); v(0)], with A0 being 2*Nt x 2.

numRows = 6*Nt;
numCols = 6;

% A = [1, dt; 0, 1]
ii_A0 = [ix_x(1), iv_x(1)];
jj_A0 = [1, 4];
vv_A0 = [1, 1];

A0 = sparse([ii_A0, ii_A0+Nt, ii_A0+2*Nt], ...
    [jj_A0, jj_A0+1, jj_A0+2], ...
    [vv_A0, vv_A0, vv_A0], numRows, numCols);



end