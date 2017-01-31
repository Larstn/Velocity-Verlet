%% Returns the Acceleration Interpolation Matrix
% Calculates and returns the matrix to interpolate the E-Field to get the
% acceleration vector at the particle positions 
function accelInterpMatrix = get_accelInterpmatrix(iix, iiy, iiz, ...
    w000, w001, w010, w011, w100, w101, w110, w111, Nx, Ny, Nz, Nt) 

    sz = [Nx Ny Nz];
    i000 = sub2ind(sz, iix, iiy, iiz);
    i001 = sub2ind(sz, iix+1, iiy, iiz);
    i010 = sub2ind(sz, iix, iiy+1, iiz);
    i011 = sub2ind(sz, iix+1, iiy+1, iiz);
    i100 = sub2ind(sz, iix, iiy, iiz+1);
    i101 = sub2ind(sz, iix+1, iiy, iiz+1);
    i110 = sub2ind(sz, iix, iiy+1, iiz+1);
    i111 = sub2ind(sz, iix+1, iiy+1, iiz+1);

    nn = 1:Nt;
    accelInterpMatrix = sparse([nn, nn, nn, nn, nn, nn, nn, nn], ...
        [i000; i001; i010; i011; i100; i101; i110; i111], ...
        [w000; w001; w010; w011; w100; w101; w110; w111], ...
        Nt, Nx*Ny*Nz);
    
    
end