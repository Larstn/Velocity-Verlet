function [result] = multipleRestriction3D(i_x , i_y, i_z, w000, w100, w010, w110, w001, w101, w011, w111, val, Nx, Ny, Nz)

result = zeros(Nx, Ny, Nz);

assert(all(w000(:) >= -1e-6))
assert(all(w100(:) >= -1e-6))
assert(all(w010(:) >= -1e-6))
assert(all(w110(:) >= -1e-6))
assert(all(w001(:) >= -1e-6))
assert(all(w101(:) >= -1e-6))
assert(all(w010(:) >= -1e-6))
assert(all(w111(:) >= -1e-6))

assert(all(i_x(:) > 0))
assert(all(i_y(:) > 0))
assert(all(i_z(:) > 0))
assert(all(i_x(:) <= Nx-1))
assert(all(i_y(:) <= Ny-1))
assert(all(i_z(:) <= Nz-1))

for ii = 1:size(i_x,1)
    
    val_ii = val(ii);
    
    result( i_x(ii)      , i_y(ii)       , i_z(ii))        = result( i_x(ii)       , i_y(ii)       , i_z(ii)) + w000(ii)*val_ii;    
    
    if(i_x(ii) < Nx)
        result( i_x(ii) + 1  , i_y(ii)       , i_z(ii))        = result( i_x(ii) + 1   , i_y(ii)       , i_z(ii)) + w100(ii)*val_ii;    end    
    if(i_y(ii) < Ny) 
        result( i_x(ii)      , i_y(ii) + 1   , i_z(ii))        = result( i_x(ii)       , i_y(ii) + 1   , i_z(ii)) + w010(ii)*val_ii;    end
    if((i_x(ii) < Nx) && (i_y(ii) < Ny))
        result( i_x(ii) + 1  , i_y(ii) + 1   , i_z(ii))        = result( i_x(ii) + 1   , i_y(ii) + 1   , i_z(ii)) + w110(ii)*val_ii;    end
    
    if(i_z(ii) < Nz)
        result( i_x(ii)      , i_y(ii)       , i_z(ii) + 1)    = result( i_x(ii)      , i_y(ii)       , i_z(ii) + 1) + w001(ii)*val_ii; end
    if((i_x(ii) < Nx) && (i_z(ii) < Nz))
        result( i_x(ii) + 1  , i_y(ii)       , i_z(ii) + 1)    = result( i_x(ii) + 1  , i_y(ii)       , i_z(ii) + 1) + w101(ii)*val_ii; end
    if((i_y(ii) < Ny) && (i_z(ii) < Nz))
        result( i_x(ii)      , i_y(ii) + 1   , i_z(ii) + 1)    = result( i_x(ii)      , i_y(ii) + 1   , i_z(ii) + 1) + w011(ii)*val_ii; end
    if((i_x(ii) < Nx) && (i_y(ii) < Ny) && (i_z(ii) < Nz))
        result( i_x(ii) + 1  , i_y(ii) + 1   , i_z(ii) + 1)    = result( i_x(ii) + 1  , i_y(ii) + 1   , i_z(ii) + 1) + w111(ii)*val_ii; end
end

end

