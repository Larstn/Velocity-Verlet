classdef Velocity_Verlet
    properties
        xyz
        resolution
        
        V
        Ex
        Ey
        Ez
        
        Nparticles
        ParticleArray
        
        Fobj            %objective function (needs to be a function handle)
        
    end 
    
    properties (SetAccess=private, GetAccess=public, Hidden)
        x_grid
        y_grid
        z_grid
        
        dx
        dy
        dz
    end
    
    properties (SetAccess=private, GetAccess=public)
        dFdEx
        dFdEy
        dFdEz
        
        dFdV
        
        Fval        
    end
    
    methods        

        function obj = Velocity_Verlet(xyz, resolution, Nparticles, ...
                ParticleArray, Fobj, varargin)
%% Inputs: Dimensions, Resolutions, Array of Particles, Function Handle Objective Function, E_fields (3) or Potential (1)            
            assert((isa(xyz,'double') && isequal(size(xyz),[3 ,2])), ...
                'Grid dimensions must be 3x2 doubles')
            obj.xyz = xyz;
            assert((isa(resolution,'double') && isequal(size(resolution)...
                ,[1 ,3])), 'Resolution must be 1x3 doubles')
            obj.resolution = resolution;
            obj.x_grid = linspace(xyz(1,1),xyz(1,2),resolution(1));                        
            obj.y_grid = linspace(xyz(2,1),xyz(2,2),resolution(2));            
            obj.z_grid = linspace(xyz(3,1),xyz(3,2),resolution(3));
            

            obj.dx = obj.x_grid(2) - obj.x_grid(1);
            obj.dy = obj.y_grid(2) - obj.y_grid(1);
            obj.dz = obj.z_grid(2) - obj.z_grid(1);
            
            obj.Nparticles = Nparticles;
            obj.ParticleArray = ParticleArray;
            obj.Fobj = Fobj;
            
            
            if nargin == 8
                
                assert((isa(varargin{1},'double') && ...
                    isequal(size(varargin{1}),[resolution(1),resolution(2),resolution(3)])), ...
                    'E_x field must be the same size as the resolution')
                assert((isa(varargin{2},'double') && ...
                    isequal(size(varargin{2}),[resolution(1),resolution(2),resolution(3)])), ...
                    'E_y field must be the same size as the resolution')
                assert((isa(varargin{3},'double')&& ...
                    isequal(size(varargin{3}),[resolution(1),resolution(2),resolution(3)])), ...
                    'E_z field must be the same size as the resolution')
                
                obj.Ex = varargin{1};
                obj.Ey = varargin{2};
                obj.Ez = varargin{3};
                
                obj.V = cumsum(varargin{1},1)*obj.dx + cumsum(varargin{2},2)*obj.dy + cumsum(varargin{3},3)*obj.dz;
                
            elseif nargin == 6
                assert((isa(varargin{1},'double') && ...
                    isequal(size(varargin{1}),[resolution(1),resolution(2),resolution(3)])), ...
                    'V field must be the same size as the resolution')
                
                obj.V   = varargin{1};
                obj.Ex  = -centeredDiff(varargin{1}, 1) / obj.dx;
                obj.Ey  = -centeredDiff(varargin{1}, 2) / obj.dy;
                obj.Ez  = -centeredDiff(varargin{1}, 3) / obj.dz;

            else
                fprintf("The input must either be three electric fields or one potential")
                assert(false)
            end 
                                    
        end
        
        function obj = calculateF_dF(obj)
            
                G_sum = 0;

                dGdEx_sum = 0*obj.Ex;
                dGdEy_sum = 0*obj.Ey;
                dGdEz_sum = 0*obj.Ez;
                nParticle = size(obj.ParticleArray,2);
                %DG_sum = zeros(6*Nt,1);
            
            for ii = 1:size(obj.ParticleArray,2)
                accelFunc = obj.accelerationFunction();
                obj.ParticleArray(ii).accelerationFunction = accelFunc(obj.Ex, obj.Ey, obj.Ez);
                obj.ParticleArray(ii) = obj.ParticleArray(ii).runTrajectory();
                
                [iix, iiy, iiz, w000, w001, w010, w011, w100, w101, w110, w111] ...
                    = trilinear_weights(obj,ii);
                    
                accelInterpMatrix = get_accelInterpmatrix(obj,iix, iiy, iiz, ...
                    w000, w001, w010, w011, w100, w101, w110, w111, ...
                    ii);
                
                [ix_x, ix_y, ix_z, ~] = obj.ParticleArray(ii).get_Index3D();
                
                [Ix, Iy, Iz] = get_I(obj, ii);
                
                [G, DG] = obj.Fobj(obj.ParticleArray(1,ii).xv);
                
                [systemMatrix, ~, accelMatrix] = ...
                    obj.velocityVerletMatrices3D(ii);


                systemMatrix = ...
                    systemMatrix*...
                    (1/((obj.ParticleArray(ii).ncharges*obj.ParticleArray(ii).e)/(obj.ParticleArray(ii).nmass*obj.ParticleArray(ii).u)));
                

                [S_p] = obj.get_PrimalS(...
                    Ix, Iy, Iz, ii, systemMatrix);

                obj.ParticleArray(1,ii).xv_dual = S_p' \ DG;
                
                [dGdEx, dGdEy, dGdEz, ~] = obj.getdGdE(ii, accelMatrix, ix_x, ix_y, ix_z, accelInterpMatrix);
                 
                G_sum = G_sum + G.^2;
        
                dGdEx_sum = dGdEx_sum + 1/nParticle*2*dGdEx*G;
                dGdEy_sum = dGdEy_sum + 1/nParticle*2*dGdEy*G;
                dGdEz_sum = dGdEz_sum + 1/nParticle*2*dGdEz*G;
                
            end
            
            obj.Fval = sqrt(1/nParticle*G_sum);
            obj.dFdEx = 0.5*dGdEx_sum ./ obj.Fval;
            obj.dFdEy = 0.5*dGdEy_sum ./ obj.Fval;
            obj.dFdEz = 0.5*dGdEz_sum ./ obj.Fval;
            
            obj.dFdV = 0*obj.dFdEx;
            
            obj.dFdV(3:end,:,:)   = obj.dFdV(3:end,:,:)     +...
                (-0.5/obj.dx)*obj.dFdEx(2:end-1,:,:);
            obj.dFdV(1:end-2,:,:) = obj.dFdV(1:end-2,:,:)   +...
                (0.5/obj.dx)*obj.dFdEx(2:end-1,:,:);

            obj.dFdV(:,3:end,:)   = obj.dFdV(:,3:end,:)     +...
                (-0.5/obj.dy)*obj.dFdEy(:,2:end-1,:);
            obj.dFdV(:,1:end-2,:) = obj.dFdV(:,1:end-2,:)   +...
                (0.5/obj.dy)*obj.dFdEy(:,2:end-1,:);

            obj.dFdV(:,:,3:end)   = obj.dFdV(:,:,3:end)     +...
                (-0.5/obj.dz)*obj.dFdEz(:,:,2:end-1);
            obj.dFdV(:,:,1:end-2) = obj.dFdV(:,:,1:end-2)   +...
                (0.5/obj.dz)*obj.dFdEz(:,:,2:end-1);
            
        end
            

        function accelFunc = accelerationFunction(obj)
%% Creation of the acceleration Function
% Creates the acceleration function with the given grid as well as the
% charge (as a multiple of the elementary charge) and the masss (as
% multiple of the electron mass). The created acceleration function uses
% the electric field in all three components as well as the position as an
% input
%%

       % elementary_charge   = 1.60217662e-19;
       % electron_mass       = 1.6605e-27;


        interpolant = @(E, x, y, z) interpn(obj.x_grid, obj.y_grid, obj.z_grid, E, x, y, z, 'linear', 0);


        accelFunc = @(Ex,Ey,Ez) @(t, xyz) [...
            interpolant(Ex, xyz(1), xyz(2), xyz(3)); ...
            interpolant(Ey, xyz(1), xyz(2), xyz(3)); ...
            interpolant(Ez, xyz(1), xyz(2), xyz(3)); ...
            ];

        end
        
        function [i_x, i_y, i_z, w000, w100, w010, w110, w001, w101, w011, w111] = trilinear_weights(obj,ii)
%% trilinear_weights    Indices and weights for trilinear interpolation
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

            %             assert(isvector(x))
            %             assert(isvector(y))
            %             assert(isvector(z))
            %             assert(isvector(x_grid))
            %             assert(isvector(y_grid))
            %             assert(isvector(z_grid))
            %             assert(isequal(size(x), size(y), size(z)))
            x = obj.ParticleArray(1,ii).xx;
            y = obj.ParticleArray(1,ii).yy;
            z = obj.ParticleArray(1,ii).zz;
                        % Coerce x_grid etc. to be row vectors if x, y, z are row vectors...
            if isrow(x) && ~isrow(obj.x_grid)
                row = @(A) reshape(A, 1, []);
                obj.x_grid = row(obj.x_grid);
                obj.y_grid = row(obj.y_grid);
                obj.z_grid = row(obj.z_grid);
            elseif iscolumn(x) && ~iscolumn(obj.x_grid)
                col = @(A) reshape(A, [], 1);
                obj.x_grid = col(obj.x_grid);
                obj.y_grid = col(obj.y_grid);
                obj.z_grid = col(obj.z_grid);
            end


            %             d_x = x_grid(2)-x_grid(1);
            %             d_y = y_grid(2)-y_grid(1);
            %             d_z = z_grid(2)-z_grid(1);
            



            log_vec = ((x >= obj.x_grid(1)) & (x <= obj.x_grid(end)) & (y >= obj.y_grid(1)) & (y <= obj.y_grid(end)) & (z >= obj.z_grid(1)) & (z <= obj.z_grid(end)));

            i_x = ones(size(x));
            i_y = ones(size(y));
            i_z = ones(size(z));

            i_x(log_vec) = floor( (x(log_vec) - obj.x_grid(1))/obj.dx) +1;
            i_y(log_vec) = floor( (y(log_vec) - obj.y_grid(1))/obj.dy) +1;
            i_z(log_vec) = floor( (z(log_vec) - obj.z_grid(1))/obj.dz) +1;

            % Handle a special case: when x == x_grid(end) we round its position DOWN
            % instead of UP, to allow all boundary values to be defined as one might
            % expect.  Now x == x_grid(1) is in the first cell AND x == x_grid(end)
            % is in the last cell.

            i_x(x == obj.x_grid(end)) = length(obj.x_grid)-1;
            i_y(y == obj.y_grid(end)) = length(obj.y_grid)-1;
            i_z(z == obj.z_grid(end)) = length(obj.z_grid)-1;

            x_c = obj.x_grid(i_x);
            y_c = obj.y_grid(i_y);
            z_c = obj.z_grid(i_z);

            % Recover weights

            wx = ( x - x_c ) ./ obj.dx;
            wy = ( y - y_c ) ./ obj.dy;
            wz = ( z - z_c ) ./ obj.dz; 

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
        
        function accelInterpMatrix = get_accelInterpmatrix(obj, iix, iiy, iiz, ...
            w000, w001, w010, w011, w100, w101, w110, w111,ii) 
%% Returns the Acceleration Interpolation Matrix
% Calculates and returns the matrix to interpolate the E-Field to get the
% acceleration vector at the particle positions  
            Nt = obj.ParticleArray(1,ii).Nt;
            
            i000 = sub2ind(obj.resolution, iix, iiy, iiz);
            i001 = sub2ind(obj.resolution, iix+1, iiy, iiz);
            i010 = sub2ind(obj.resolution, iix, iiy+1, iiz);
            i011 = sub2ind(obj.resolution, iix+1, iiy+1, iiz);
            i100 = sub2ind(obj.resolution, iix, iiy, iiz+1);
            i101 = sub2ind(obj.resolution, iix+1, iiy, iiz+1);
            i110 = sub2ind(obj.resolution, iix, iiy+1, iiz+1);
            i111 = sub2ind(obj.resolution, iix+1, iiy+1, iiz+1);

            nn = 1:Nt;
            accelInterpMatrix = sparse([nn, nn, nn, nn, nn, nn, nn, nn], ...
                [i000; i001; i010; i011; i100; i101; i110; i111], ...
                [w000; w001; w010; w011; w100; w101; w110; w111], ...
                Nt, prod(obj.resolution));

        end
        

        function [Ix, Iy, Iz] = get_I(obj, ii)
%% Calculate Ix, Iy, Iz (Derivative of Trilinear Weights)
% Calculates the derivative of the trilinear weights with respect to x,y,z

            %xv = obj.ParticleArray(1,ii).xv;
            %Nt = obj.ParticleArray(1,ii).Nt;
            
            Ix = obj.trilinear_weights_I_x(ii);
            Iy = obj.trilinear_weights_I_y(ii);
            Iz = obj.trilinear_weights_I_z(ii);

        end 
        function [partial_x_accelInterpMatrix, i_x, i_y, i_z, w000, w001, w010, w011, w100, w101, w110, w111] = trilinear_weights_I_x(obj,ii)
%% trilinear_weights    Indices and weights for trilinear interpolation
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


            x = obj.ParticleArray(1,ii).xx;
            y = obj.ParticleArray(1,ii).yy;
            z = obj.ParticleArray(1,ii).zz;
            Nt = obj.ParticleArray(ii).Nt;
            assert(isvector(x))
            assert(isvector(y))
            assert(isvector(z))
            assert(isvector(obj.x_grid))
            assert(isvector(obj.y_grid))
            assert(isvector(obj.z_grid))
            assert(isequal(size(x), size(y), size(z)))
            % Coerce x_grid etc. to be row vectors if x, y, z are row vectors...
            if isrow(x) && ~isrow(obj.x_grid)
                row = @(A) reshape(A, 1, []);
                obj.x_grid = row(obj.x_grid);
                obj.y_grid = row(obj.y_grid);
                obj.z_grid = row(obj.z_grid);
            elseif iscolumn(x) && ~iscolumn(obj.x_grid)
                col = @(A) reshape(A, [], 1);
                obj.x_grid = col(obj.x_grid);
                obj.y_grid = col(obj.y_grid);
                obj.z_grid = col(obj.z_grid);
            end


            

            % TODO: what does "log" mean in "log_vec"? - It was supposed to mean
            % "Logic" just because it is a binary vector
            log_vec = ((x >= obj.x_grid(1)) & (x <= obj.x_grid(end)) & (y >= obj.y_grid(1)) & (y <= obj.y_grid(end)) & (z >= obj.z_grid(1)) & (z <= obj.z_grid(end)));

            i_x = ones(size(x));
            i_y = ones(size(y));
            i_z = ones(size(z));

            i_x(log_vec) = floor( (x(log_vec) - obj.x_grid(1))/obj.dx) +1;
            i_y(log_vec) = floor( (y(log_vec) - obj.y_grid(1))/obj.dy) +1;
            i_z(log_vec) = floor( (z(log_vec) - obj.z_grid(1))/obj.dz) +1;

            % Handle a special case: when x == obj.x_grid(end) we round its position DOWN
            % instead of UP, to allow all boundary values to be defined as one might
            % expect.  Now x == obj.x_grid(1) is in the first cell AND x == obj.x_grid(end)
            % is in the last cell.

            i_x(x == obj.x_grid(end)) = length(obj.x_grid)-1;
            i_y(y == obj.y_grid(end)) = length(obj.y_grid)-1;
            i_z(z == obj.z_grid(end)) = length(obj.z_grid)-1;

            x_c = obj.x_grid(i_x);
            y_c = obj.y_grid(i_y);
            z_c = obj.z_grid(i_z);

            % Recover weights

            wx = ( 1 ) ./ obj.dx;
            wy = ( y - y_c ) ./ obj.dy;
            wz = ( z - z_c ) ./ obj.dz; 

            w000 = (-wx).*(1-wy).*(1-wz);
            w100 = wx.*(1-wy).*(1-wz);
            w010 = (-wx).*wy.*(1-wz);
            w110 = wx.*wy.*(1-wz);
            w001 = (-wx).*(1-wy).*wz;
            w101 = wx.*(1-wy).*wz;
            w011 = (-wx).*wy.*wz;
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



           
            i000 = sub2ind(obj.resolution, i_x, i_y, i_z);
            i001 = sub2ind(obj.resolution, i_x, i_y, i_z+1);
            i010 = sub2ind(obj.resolution, i_x, i_y+1, i_z);
            i011 = sub2ind(obj.resolution, i_x, i_y+1, i_z+1);
            i100 = sub2ind(obj.resolution, i_x+1, i_y, i_z);
            i101 = sub2ind(obj.resolution, i_x+1, i_y, i_z+1);
            i110 = sub2ind(obj.resolution, i_x+1, i_y+1, i_z);
            i111 = sub2ind(obj.resolution, i_x+1, i_y+1, i_z+1);


            nn = 1:Nt;
            partial_x_accelInterpMatrix = sparse([nn, nn, nn, nn, nn, nn, nn, nn], ...
                [i000; i001; i010; i011; i100; i101; i110; i111], ...
                [w000; w001; w010; w011; w100; w101; w110; w111], ...
                Nt, prod(obj.resolution));

        end     
        

function [partial_x_accelInterpMatrix, i_x, i_y, i_z, w000, w001, w010, w011, w100, w101, w110, w111] = trilinear_weights_I_y(obj, ii)
%% trilinear_weights    Indices and weights for trilinear interpolation
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
            

            x = obj.ParticleArray(1,ii).xx;
            y = obj.ParticleArray(1,ii).yy;
            z = obj.ParticleArray(1,ii).zz;
            Nt = obj.ParticleArray(ii).Nt;
            
            assert(isvector(x))
            assert(isvector(y))
            assert(isvector(z))
            assert(isvector(obj.x_grid))
            assert(isvector(obj.y_grid))
            assert(isvector(obj.z_grid))
            assert(isequal(size(x), size(y), size(z)))

            % Coerce obj.x_grid etc. to be row vectors if x, y, z are row vectors...
            if isrow(x) && ~isrow(obj.x_grid)
                row = @(A) reshape(A, 1, []);
                obj.x_grid = row(obj.x_grid);
                obj.y_grid = row(obj.y_grid);
                obj.z_grid = row(obj.z_grid);
            elseif iscolumn(x) && ~iscolumn(obj.x_grid)
                col = @(A) reshape(A, [], 1);
                obj.x_grid = col(obj.x_grid);
                obj.y_grid = col(obj.y_grid);
                obj.z_grid = col(obj.z_grid);
            end


            % TODO: what does "log" mean in "log_vec"? - It was supposed to mean
            % "Logic" just because it is a binary vector
            log_vec = ((x >= obj.x_grid(1)) & (x <= obj.x_grid(end)) & (y >= obj.y_grid(1)) & (y <= obj.y_grid(end)) & (z >= obj.z_grid(1)) & (z <= obj.z_grid(end)));

            i_x = ones(size(x));
            i_y = ones(size(y));
            i_z = ones(size(z));

            i_x(log_vec) = floor( (x(log_vec) - obj.x_grid(1))/obj.dx) +1;
            i_y(log_vec) = floor( (y(log_vec) - obj.y_grid(1))/obj.dy) +1;
            i_z(log_vec) = floor( (z(log_vec) - obj.z_grid(1))/obj.dz) +1;

            % Handle a special case: when x == obj.x_grid(end) we round its position DOWN
            % instead of UP, to allow all boundary values to be defined as one might
            % expect.  Now x == obj.x_grid(1) is in the first cell AND x == obj.x_grid(end)
            % is in the last cell.

            i_x(x == obj.x_grid(end)) = length(obj.x_grid)-1;
            i_y(y == obj.y_grid(end)) = length(obj.y_grid)-1;
            i_z(z == obj.z_grid(end)) = length(obj.z_grid)-1;

            x_c = obj.x_grid(i_x);
            y_c = obj.y_grid(i_y);
            z_c = obj.z_grid(i_z);

            % Recover weights

            wx = ( x - x_c ) ./ obj.dx;
            wy = ( 1 ) ./ obj.dy;
            wz = ( z - z_c ) ./ obj.dz; 

            w000 = (1-wx).*(-wy).*(1-wz);
            w100 = wx.*(-wy).*(1-wz);
            w010 = (1-wx).*wy.*(1-wz);
            w110 = wx.*wy.*(1-wz);
            w001 = (1-wx).*(-wy).*wz;
            w101 = wx.*(-wy).*wz;
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



      
            i000 = sub2ind(obj.resolution, i_x, i_y, i_z);
            i001 = sub2ind(obj.resolution, i_x, i_y, i_z+1);
            i010 = sub2ind(obj.resolution, i_x, i_y+1, i_z);
            i011 = sub2ind(obj.resolution, i_x, i_y+1, i_z+1);
            i100 = sub2ind(obj.resolution, i_x+1, i_y, i_z);
            i101 = sub2ind(obj.resolution, i_x+1, i_y, i_z+1);
            i110 = sub2ind(obj.resolution, i_x+1, i_y+1, i_z);
            i111 = sub2ind(obj.resolution, i_x+1, i_y+1, i_z+1);


            % This Part does not work for Nt = 2
            nn = 1:Nt;
            partial_x_accelInterpMatrix = sparse([nn, nn, nn, nn, nn, nn, nn, nn], ...
                [i000; i001; i010; i011; i100; i101; i110; i111], ...
                [w000; w001; w010; w011; w100; w101; w110; w111], ...
                Nt, prod(obj.resolution));

            % checkClose = @(a, b) assert(norm(a-b) < 1e-6);

            % 
            % A = full(partial_x_accelInterpMatrix);
            % checkClose([w000(1), w001(1), w010(1), w011(1), w100(1), w101(1), w110(1), w111(1)], [A(1,1), A(2,2), A(3,12), A(1, 13), A(2,122 ), A(3,123 ), A(1, 133), A(2, 134)]);
            % 
            % checkClose([w000(2), w001(2), w010(2), w011(2), w100(2), w101(2), w110(2), w111(2)], [A(3,1), A(1,2), A(2,12), A(3, 13), A(1,122 ), A(2,123 ), A(3, 133), A(1, 134)]);
            % 
            % checkClose([w000(3), w001(3), w010(3), w011(3), w100(3), w101(3), w110(3), w111(3)], [A(2,1), A(3,2), A(1,12), A(2, 13), A(3,122 ), A(1,123 ), A(2, 133), A(3, 134)]);
            % 

        end
        

        function [partial_x_accelInterpMatrix, i_x, i_y, i_z, w000, w001, w010, w011, w100, w101, w110, w111] = trilinear_weights_I_z(obj,ii)
%% trilinear_weights    Indices and weights for trilinear interpolation
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
            x = obj.ParticleArray(1,ii).xx;
            y = obj.ParticleArray(1,ii).yy;
            z = obj.ParticleArray(1,ii).zz;
            Nt = obj.ParticleArray(ii).Nt;
            
            assert(isvector(x))
            assert(isvector(y))
            assert(isvector(z))
            assert(isvector(obj.x_grid))
            assert(isvector(obj.y_grid))
            assert(isvector(obj.z_grid))
            assert(isequal(size(x), size(y), size(z)))

            % Coerce obj.x_grid etc. to be row vectors if x, y, z are row vectors...
            if isrow(x) && ~isrow(obj.x_grid)
                row = @(A) reshape(A, 1, []);
                obj.x_grid = row(obj.x_grid);
                obj.y_grid = row(obj.y_grid);
                obj.z_grid = row(obj.z_grid);
            elseif iscolumn(x) && ~iscolumn(obj.x_grid)
                col = @(A) reshape(A, [], 1);
                obj.x_grid = col(obj.x_grid);
                obj.y_grid = col(obj.y_grid);
                obj.z_grid = col(obj.z_grid);
            end



            % TODO: what does "log" mean in "log_vec"? - It was supposed to mean
            % "Logic" just because it is a binary vector
            log_vec = ((x >= obj.x_grid(1)) & (x <= obj.x_grid(end)) & (y >= obj.y_grid(1)) & (y <= obj.y_grid(end)) & (z >= obj.z_grid(1)) & (z <= obj.z_grid(end)));

            i_x = ones(size(x));
            i_y = ones(size(y));
            i_z = ones(size(z));

            i_x(log_vec) = floor( (x(log_vec) - obj.x_grid(1))/obj.dx) +1;
            i_y(log_vec) = floor( (y(log_vec) - obj.y_grid(1))/obj.dy) +1;
            i_z(log_vec) = floor( (z(log_vec) - obj.z_grid(1))/obj.dz) +1;

            % Handle a special case: when x == obj.x_grid(end) we round its position DOWN
            % instead of UP, to allow all boundary values to be defined as one might
            % expect.  Now x == obj.x_grid(1) is in the first cell AND x == obj.x_grid(end)
            % is in the last cell.

            i_x(x == obj.x_grid(end)) = length(obj.x_grid)-1;
            i_y(y == obj.y_grid(end)) = length(obj.y_grid)-1;
            i_z(z == obj.z_grid(end)) = length(obj.z_grid)-1;

            x_c = obj.x_grid(i_x);
            y_c = obj.y_grid(i_y);
            z_c = obj.z_grid(i_z);

            % Recover weights

            wx = ( x - x_c ) ./ obj.dx;
            wy = ( y - y_c ) ./ obj.dy;
            wz = ( 1 ) ./ obj.dz; 

            w000 = (1-wx).*(1-wy).*(-wz);
            w100 = wx.*(1-wy).*(-wz);
            w010 = (1-wx).*wy.*(-wz);
            w110 = wx.*wy.*(-wz);
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

            i000 = sub2ind(obj.resolution, i_x, i_y, i_z);
            i001 = sub2ind(obj.resolution, i_x, i_y, i_z+1);
            i010 = sub2ind(obj.resolution, i_x, i_y+1, i_z);
            i011 = sub2ind(obj.resolution, i_x, i_y+1, i_z+1);
            i100 = sub2ind(obj.resolution, i_x+1, i_y, i_z);
            i101 = sub2ind(obj.resolution, i_x+1, i_y, i_z+1);
            i110 = sub2ind(obj.resolution, i_x+1, i_y+1, i_z);
            i111 = sub2ind(obj.resolution, i_x+1, i_y+1, i_z+1);


            % This Part does not work for Nt = 2
            nn = 1:Nt;
            partial_x_accelInterpMatrix = sparse([nn, nn, nn, nn, nn, nn, nn, nn], ...
                [i000; i001; i010; i011; i100; i101; i110; i111], ...
                [w000; w001; w010; w011; w100; w101; w110; w111], ...
                Nt, prod(obj.resolution));

            % checkClose = @(a, b) assert(norm(a-b) < 1e-6);

            % 
            % A = full(partial_x_accelInterpMatrix);
            % checkClose([w000(1), w001(1), w010(1), w011(1), w100(1), w101(1), w110(1), w111(1)], [A(1,1), A(2,2), A(3,12), A(1, 13), A(2,122 ), A(3,123 ), A(1, 133), A(2, 134)]);
            % 
            % checkClose([w000(2), w001(2), w010(2), w011(2), w100(2), w101(2), w110(2), w111(2)], [A(3,1), A(1,2), A(2,12), A(3, 13), A(1,122 ), A(2,123 ), A(3, 133), A(1, 134)]);
            % 
            % checkClose([w000(3), w001(3), w010(3), w011(3), w100(3), w101(3), w110(3), w111(3)], [A(2,1), A(3,2), A(1,12), A(2, 13), A(3,122 ), A(1,123 ), A(2, 133), A(3, 134)]);
            % 

        end

        function [systemMatrix, A0, B] = velocityVerletMatrices3D(obj,ii)
%% velocityVerletMatrices3D    Matrices and indexing for 3D velocity Verlet method
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

            %Nt = length(ts);
            
            %dt = (t_end - t_start)./(Nt-1);

            [ix_x, ix_y, ix_z, iv_x, iv_y, iv_z] = get_Index3D(obj.ParticleArray(ii).Nt);

            ia_x = ix_x;
            ia_y = ix_y;
            ia_z = ix_z;


            unos = @(n) ones(1,n);

            % ------- Diagonal
    

            ii_diag = [ix_x(1:obj.ParticleArray(ii).Nt), iv_x(1:obj.ParticleArray(ii).Nt)];
            jj_diag = [ix_x(1:obj.ParticleArray(ii).Nt), iv_x(1:obj.ParticleArray(ii).Nt)];
            vv_diag = [-unos(obj.ParticleArray(ii).Nt), -unos(obj.ParticleArray(ii).Nt)];



            ii_A = [ix_x(2:obj.ParticleArray(ii).Nt), ix_x(2:obj.ParticleArray(ii).Nt), iv_x(2:obj.ParticleArray(ii).Nt)];
            jj_A = [ix_x(1:obj.ParticleArray(ii).Nt-1), iv_x(1:obj.ParticleArray(ii).Nt-1), iv_x(1:obj.ParticleArray(ii).Nt-1)];
            vv_A = [unos(obj.ParticleArray(ii).Nt-1), obj.ParticleArray(ii).dt*unos(obj.ParticleArray(ii).Nt-1), unos(obj.ParticleArray(ii).Nt-1)];



            numRows = 6*obj.ParticleArray(ii).Nt;
            numCols = 6*obj.ParticleArray(ii).Nt;

            systemMatrix = sparse(...
                [ii_A, ii_diag, ii_A+obj.ParticleArray(ii).Nt, ii_diag+obj.ParticleArray(ii).Nt, ii_A+2*obj.ParticleArray(ii).Nt, ii_diag+2*obj.ParticleArray(ii).Nt], ...
                [jj_A, jj_diag, jj_A+obj.ParticleArray(ii).Nt, jj_diag+obj.ParticleArray(ii).Nt, jj_A+2*obj.ParticleArray(ii).Nt, jj_diag+2*obj.ParticleArray(ii).Nt], ...
                [vv_A, vv_diag, vv_A, vv_diag, vv_A, vv_diag],...
                numRows, numCols);


            % Right-hand side


            ii_B = [ix_x(2:obj.ParticleArray(ii).Nt), iv_x(2:obj.ParticleArray(ii).Nt), iv_x(2:obj.ParticleArray(ii).Nt)];
            jj_B = [ia_x(1:obj.ParticleArray(ii).Nt-1), ia_x(1:obj.ParticleArray(ii).Nt-1), ia_x(2:obj.ParticleArray(ii).Nt)];
            vv_B = [0.5*obj.ParticleArray(ii).dt*obj.ParticleArray(ii).dt*unos(obj.ParticleArray(ii).Nt-1), 0.5*obj.ParticleArray(ii).dt*unos(obj.ParticleArray(ii).Nt-1), 0.5*obj.ParticleArray(ii).dt*unos(obj.ParticleArray(ii).Nt-1)];

            numRows = 6*obj.ParticleArray(ii).Nt;
            numCols = 3*obj.ParticleArray(ii).Nt;
            B = sparse([ii_B, ii_B+obj.ParticleArray(ii).Nt, ii_B+2*obj.ParticleArray(ii).Nt],...
                [jj_B, jj_B+obj.ParticleArray(ii).Nt, jj_B+2*obj.ParticleArray(ii).Nt], ...
                [vv_B, vv_B, vv_B], numRows, numCols);
            %figure(3); clf
            %imagesc(B)

            % ------- Initial conditions
            % The first timestep needs to take initial conditions:
            % [x(1); v(1)] = A*[x(0); v(0)] + B*[a(0); a(1)]
            % So we need an RHS vector of the right size, 2*Nt x 1.
            % Let's get RHS = A0 * [x(0); v(0)], with A0 being 2*Nt x 2.

            numRows = 6*obj.ParticleArray(ii).Nt;
            numCols = 6;

            % A = [1, dt; 0, 1]
            ii_A0 = [ix_x(1), iv_x(1)];
            jj_A0 = [1, 4];
            vv_A0 = [1, 1];

            A0 = sparse([ii_A0, ii_A0+obj.ParticleArray(ii).Nt, ii_A0+2*obj.ParticleArray(ii).Nt], ...
                [jj_A0, jj_A0+1, jj_A0+2], ...
                [vv_A0, vv_A0, vv_A0], numRows, numCols);



        end
        

        function [S_p] = get_PrimalS(obj, Ix, Iy, Iz, ii, systemMatrix)
%% Calculation of the Primal Matrix
% This function calculates the matrix of the Primal System (for a
% derivation with respect to the electric field).
                D_I_1          = [Ix*obj.Ex(:); Ix*obj.Ey(:); Ix*obj.Ez(:)]; % [kg*m*s^-2*A^-1]
                D_I_2          = [Iy*obj.Ex(:); Iy*obj.Ey(:); Iy*obj.Ez(:)];
                D_I_3          = [Iz*obj.Ex(:); Iz*obj.Ey(:); Iz*obj.Ez(:)];

                grad_E_1 = zeros(1,3*obj.ParticleArray(ii).Nt);
                grad_E_2 = zeros(1,3*obj.ParticleArray(ii).Nt);
                grad_E_3 = zeros(1,3*obj.ParticleArray(ii).Nt);

                Bp = B_matrix_primal(obj.ParticleArray(ii).t_vec,grad_E_1, grad_E_2, grad_E_3, D_I_1, D_I_2, D_I_3);


                S_p = systemMatrix + Bp;
        end


        function [dGdEx, dGdEy, dGdEz, dGda] = getdGdE(obj, ii, accelMatrix, ...
            ix_x, ix_y, ix_z, accelInterpMatrix)
%% dGdEx, dGdEy, dGdEz Calculation
% This function calculates the derivatives of the objective function
% towards the components of the electric field at each grid point

            xv_dual = obj.ParticleArray(ii).xv_dual;
            dGda = xv_dual' * (-accelMatrix);

            dGdax = dGda(ix_x);
            dGdEx = reshape(dGdax * accelInterpMatrix, obj.resolution);

            dGday = dGda(ix_y);
            dGdEy = reshape(dGday * accelInterpMatrix, obj.resolution);

            dGdaz = dGda(ix_z);
            dGdEz = reshape(dGdaz * accelInterpMatrix, obj.resolution);


        end
        
        function plotdFdEx(obj)
            
            figure(1010)
            imagesc(obj.x_grid,obj.y_grid, obj.dFdEx(:,:,1)')
            axis xy image
            colorbar
            xlabel('x [m]')
            ylabel('y [m]')
        end
    end
    
end
