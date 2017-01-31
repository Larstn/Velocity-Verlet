function [xv, accelerations] = velocityVerlet3D(ts, xv0, accelFunc)
% velocityVerlet3D     Calculate trajectory with driving force
%
% [xv, accel, ix_x, ix_y, ix_z, iv_x, iv_y, iv_z] = velocityVerlet3D(ts, xv0, accelFunc)
%
% accelFunc(t, x) returns the particle acceleration at one time and
% position, as a 3-vector [ax; ay; az].
%
% xv0 = [x0; y0; z0; vx0; vy0; vz0] are the initial conditions.
% ts should be uniformly-spaced but this is not checked.
%
% xv is the position and velocity at times ts(2:end).
% accel is the acceleration applied at times ts(1:end).
% Note that these are not the same length.
%
% accel should not be needed in practice but may be useful for debugging.
% (PCH uses it to compare values from accelFunc with trilinear
% interpolation done "by hand".)

Nt = length(ts);
dt = ts(2) - ts(1);


% System state: position and velocity
xs = zeros(Nt, 3);
vs = zeros(Nt, 3);
xs(1,1:3) = xv0(1:3); % m/s
vs(1,1:3) = xv0(4:6); % m/s

accelerations = zeros(Nt, 3);

currentAcceleration = accelFunc(ts(1), xs(1,1:3));
assert(isequal(size(currentAcceleration), [3, 1]));

accelerations(1,1:3) = currentAcceleration;

% Integrate!
for nn = 1:Nt-1
    
    xs(nn+1,:) = xs(nn,:) + vs(nn,:)*dt + 0.5*dt*dt*currentAcceleration';
    nextAcceleration = accelFunc(ts(nn+1), xs(nn+1,:));
    vs(nn+1,:) = vs(nn,:) + 0.5*dt*(currentAcceleration + nextAcceleration)';
    accelerations(nn+1,1:3) = nextAcceleration;
    
    currentAcceleration = nextAcceleration;
end

xv = [xs(:,1); xs(:,2); xs(:,3);
    vs(:,1); vs(:,2); vs(:,3)];



