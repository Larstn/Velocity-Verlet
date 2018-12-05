classdef cylindricalVelocity_VerletArray
    properties
        xyz
        resolution
        
        V
        Ex
        Er
        Ez
        
        Nparticles
        ParticleArray
        
        Fobj            %objective function (needs to be a function handle)
        accelFunc
    end 
    
    properties (SetAccess=private, GetAccess=public, Hidden)
        x_grid
        y_grid
        r_grid
        z_grid
        
        dx
        dy
        dr
        dz
        
        log_vec
        i_x
        i_y
        i_z
        x_c
        y_c
        z_c
    end
    
    properties (SetAccess=private, GetAccess=public)
        dFdEx
        dFdEr
        dFdEz
        
        dFdV
        
        Fval
        
        
    end
    
    methods        

        function obj = cylindricalVelocity_VerletArray(xyz, resolution, ParticleArray, Fobj,varargin)
%% Inputs: Dimensions, Resolutions, Array of Particles, Function Handle Objective Function, E_fields (3) or Potential (1)
            assert((isa(xyz,'double') && isequal(size(xyz),[3 ,2])), ...
                'Grid dimensions must be 3x2 doubles')
            obj.xyz = xyz;
            assert((isa(resolution,'double') && isequal(size(resolution)...
                ,[1 ,3])), 'Resolution must be 1x3 doubles')
            %obj.resolution = resolution;
            obj.x_grid = linspace(xyz(1,1),xyz(1,2),resolution(1));                        
            obj.y_grid = linspace(xyz(2,1),xyz(2,2),resolution(2));            
            obj.z_grid = linspace(xyz(3,1),xyz(3,2),resolution(3));
            obj.r_grid = [-obj.y_grid(end:-1:2) obj.y_grid];
            obj.resolution = [length(obj.x_grid), length(obj.r_grid), length(obj.z_grid)];
            
            
            obj.dx = obj.x_grid(2) - obj.x_grid(1);
            obj.dy = obj.y_grid(2) - obj.y_grid(1);
            obj.dz = obj.z_grid(2) - obj.z_grid(1);
            obj.dr = obj.dy;
            
            obj.Nparticles = length(ParticleArray);
            obj.ParticleArray = ParticleArray;
            obj.Fobj = Fobj;
            
            
            if nargin == 7
                
                assert((isa(varargin{1},'double') && ...
                    isequal(size(varargin{1}),[resolution(1),resolution(2)])), ...
                    'E_x field must be the same size as the resolution')
                assert((isa(varargin{2},'double') && ...
                    isequal(size(varargin{2}),[resolution(1),resolution(2)])), ...
                    'E_y field must be the same size as the resolution')
                assert((isa(varargin{3},'double')&& ...
                    isequal(size(varargin{3}),[resolution(1),resolution(2)])), ...
                    'E_z field must be the same size as the resolution')
                
                E_x     = cat(3,varargin{1},varargin{1});
                obj.Ex = [E_x(:,end:-1:2,:) E_x(:,:,:)];
                E_y     = cat(3,varargin{2},varargin{2});
                obj.Er = [E_y(:,end:-1:2,:) E_y(:,:,:)];
                E_z     = cat(3,varargin{2},varargin{2});
                obj.Ez = [E_z(:,end:-1:2,:) E_z(:,:,:)];
                
                obj.V = cumsum(varargin{1},1)*obj.dx + cumsum(varargin{2},2)*obj.dy + cumsum(varargin{3},3)*obj.dz;
                
            elseif nargin == 5
                assert((isa(varargin{1},'double') && ...
                    isequal(size(varargin{1}),[resolution(1),resolution(2)])), ...
                    'V field must be the same size as the resolution')
                U       = cat(3,varargin{1},varargin{1});
                obj.V   = [U(:,end:-1:2,:) U(:,:,:)];
                obj.Ex  = -centeredDiff(obj.V, 1) / obj.dx;
                obj.Er  = -centeredDiff(obj.V, 2) / obj.dr;
                obj.Ez  = -centeredDiff(obj.V, 3) / obj.dz;

            else
                fprintf("The input must either be three electric fields or one potential")
                assert(false)
            end 
                                    
        end
        
        function obj = calculateF(obj)
            obj.accelFunc = obj.accelerationFunction();
            obj.ParticleArray.accelerationFunction = obj.accelFunc(obj.Ex, obj.Er, obj.Ez);
            obj.ParticleArray = obj.ParticleArray.runTrajectory();
        end
        
        function obj = calculateF_dF(obj)
            
            tic
            fprintf('Forward VV')
            obj = obj.calculateF();
            toc
            fprintf('\n')
            fprintf('Dual VV')
            tic

            obj = obj.set_logvec;

           % [w000, w001, w010, w011, w100, w101, w110, w111] ...
           %     = trilinear_weights(obj);

            %[accelInterpMatrix] = obj.get_accelInterpmatrix(w000, w001, w010, w011, w100, w101, w110, w111);

            [D_I_1, D_I_2, D_I_3] = get_I(obj);

            [G_sum, DG] = obj.Fobj(obj.ParticleArray.xv);

            DG = reshape(DG, [], 1);

            % Matricies only depent on the time vector and thus, they
            % are the same for every particle
            [systemMatrix, ~, accelMatrix] = ...
                obj.velocityVerletMatrices3D;

            systemMatrix = ...
                systemMatrix.*...
                (1/((obj.ParticleArray.ncharges*obj.ParticleArray.e)/(obj.ParticleArray.nmass*obj.ParticleArray.u)));

            [obj, S_p] = obj.get_PrimalS2(...
                D_I_1, D_I_2, D_I_3, systemMatrix);

            xv_dual_vec = S_p' \ DG;
            
            obj.ParticleArray.xv_dual = reshape(xv_dual_vec , [], obj.ParticleArray.N_particles);

            [dGdEi_cell] = obj.getdGdE(accelMatrix, xv_dual_vec);

            obj.Fval = G_sum;            
            
            dGdV = 0*dGdEi_cell{1};
            
            dGdV(3:end,:,:)   = dGdV(3:end,:,:)     +...
                (-0.5/obj.dx)*dGdEi_cell{1}(2:end-1,:,:);
            dGdV(1:end-2,:,:) = dGdV(1:end-2,:,:)   +...
                (0.5/obj.dx)*dGdEi_cell{1}(2:end-1,:,:);

            dGdV(:,3:end,:)   = dGdV(:,3:end,:)     +...
                (-0.5/obj.dy)*dGdEi_cell{2}(:,2:end-1,:);
            dGdV(:,1:end-2,:) = dGdV(:,1:end-2,:)   +...
                (0.5/obj.dy)*dGdEi_cell{2}(:,2:end-1,:);

            dGdV(:,:,3:end)   = dGdV(:,:,3:end)     +...
                (-0.5/obj.dz)*dGdEi_cell{3}(:,:,2:end-1);
            dGdV(:,:,1:end-2) = dGdV(:,:,1:end-2)   +...
                (0.5/obj.dz)*dGdEi_cell{3}(:,:,2:end-1);           
           
            dFdEx_full = sum(dGdEi_cell{1},3);
            dFdEr_full = sum(dGdEi_cell{2},3);
            dFdEz_full = sum(dGdEi_cell{3},3);
            
            dFdV_full = sum(dGdV,3);
            
            obj.dFdEx = obj.reflect_back(dFdEx_full);
            obj.dFdEr = obj.reflect_back(dFdEr_full);
            obj.dFdEz = obj.reflect_back(dFdEz_full);
           
            obj.dFdV = obj.reflect_back(dFdV_full);
            toc
            fprintf('\n')
        end
            

        function accelFunc = accelerationFunction(obj)
%% Creation of the acceleration Function
% Creates the acceleration function with the given grid as well as the
% charge (as a multiple of the elementary charge) and the masss (as
% multiple of the electron mass). The created acceleration function uses
% the electric field in all three components as well as the position as an
% input
%%



        interpolant = @(E, x, y, z) interpn(obj.x_grid, obj.r_grid, obj.z_grid, E, x, y, z, 'linear', 0);
        %interpolant = @(E, x, y, z) interplars(obj.x_grid, obj.r_grid, obj.z_grid, E, x, y, z);

        accelFunc = @(Ex,Ey,Ez) @(t, xyz) [...
            interpolant(Ex, xyz(1,:), xyz(2,:), xyz(3,:)); ...
            interpolant(Ey, xyz(1,:), xyz(2,:), xyz(3,:)); ...
            interpolant(Ez, xyz(1,:), xyz(2,:), xyz(3,:)); ...
            ];

        end
        
        
        
        
%         function accelFunc = accelerationFunction2(obj)
% %% Creation of the acceleration Function
% % Creates the acceleration function with the given grid as well as the
% % charge (as a multiple of the elementary charge) and the masss (as
% % multiple of the electron mass). The created acceleration function uses
% % the electric field in all three components as well as the position as an
% % input
% %%
% 
% 
% 
%         interpolant = @(E, x, y, z) interpn(obj.x_grid, obj.r_grid, obj.z_grid, E, x, y, z, 'linear', 0);
%         %interpolant = @(E, x, y, z) interplars(obj.x_grid, obj.r_grid, obj.z_grid, E, x, y, z);
% 
%         accelFunc = @(Ex,Ey,Ez) @(t, xyz) [...
%             interpolant(Ex, xyz(1), xyz(2), xyz(3)); ...
%             interpolant(Ey, xyz(1), xyz(2), xyz(3)); ...
%             interpolant(Ez, xyz(1), xyz(2), xyz(3)); ...
%             ];

%         end
        
        function [w000, w100, w010, w110, w001, w101, w011, w111] = trilinear_weights(obj)
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

         
            x = obj.ParticleArray.xx;
            y = obj.ParticleArray.yy;
            z = obj.ParticleArray.zz;


            % Recover weights
            
            wx = ( x - obj.x_c ) ./ obj.dx;
            wy = ( y - obj.y_c ) ./ obj.dy;
            wz = ( z - obj.z_c ) ./ obj.dz; 

            w000 = (1-wx).*(1-wy).*(1-wz);
            w100 = wx.*(1-wy).*(1-wz);
            w010 = (1-wx).*wy.*(1-wz);
            w110 = wx.*wy.*(1-wz);
            w001 = (1-wx).*(1-wy).*wz;
            w101 = wx.*(1-wy).*wz;
            w011 = (1-wx).*wy.*wz;
            w111 = wx.*wy.*wz;

            w000(~obj.log_vec) = 0;
            w100(~obj.log_vec) = 0;
            w010(~obj.log_vec) = 0;
            w110(~obj.log_vec) = 0;

            w001(~obj.log_vec) = 0;
            w101(~obj.log_vec) = 0;
            w011(~obj.log_vec) = 0;
            w111(~obj.log_vec) = 0;

            assert(isequal(size(obj.i_x), size(x), size(obj.i_y), size(y), size(obj.i_z), size(z), ...
                size(w000), size(w001), size(w010), size(w011), ...
                size(w100), size(w101), size(w110), size(w111)));



        end
        
        function [accelInterpMatrix] = get_accelInterpmatrix(...
            obj, w000, w001, w010, w011, w100, w101, w110, w111) 
%% Returns the Acceleration Interpolation Matrix
% Calculates and returns the matrix to interpolate the E-Field to get the
% acceleration vector at the particle positions  
            Nt = obj.ParticleArray.Nt;
            
            i000 = sub2ind(obj.resolution, obj.i_x, obj.i_y, obj.i_z);
            i001 = sub2ind(obj.resolution, obj.i_x+1, obj.i_y, obj.i_z);
            i010 = sub2ind(obj.resolution, obj.i_x, obj.i_y+1, obj.i_z);
            i011 = sub2ind(obj.resolution, obj.i_x+1, obj.i_y+1, obj.i_z);
            i100 = sub2ind(obj.resolution, obj.i_x, obj.i_y, obj.i_z+1);
            i101 = sub2ind(obj.resolution, obj.i_x+1, obj.i_y, obj.i_z+1);
            i110 = sub2ind(obj.resolution, obj.i_x, obj.i_y+1, obj.i_z+1);
            i111 = sub2ind(obj.resolution, obj.i_x+1, obj.i_y+1, obj.i_z+1);

            nn = 1:Nt;
            idx1_vec = [];
            idx2_vec = [];
            val_vec = [];
            
            for ii = 1:obj.ParticleArray.N_particles
                
               idx1 = [nn, nn, nn, nn, nn, nn, nn, nn] + Nt*(ii-1);
               idx2 = [i000(:,ii); i001(:,ii); i010(:,ii); i011(:,ii); i100(:,ii); i101(:,ii); i110(:,ii); i111(:,ii)] + prod(obj.resolution)*(ii-1);
               val  = [w000(:,ii); w001(:,ii); w010(:,ii); w011(:,ii); w100(:,ii); w101(:,ii); w110(:,ii); w111(:,ii)];
               
               idx1_vec = [idx1_vec, idx1(:)'];
               idx2_vec = [idx2_vec, idx2(:)'];
               val_vec  = [val_vec, val(:)'];
            
            end
            
            accelInterpMatrix = sparse(idx1_vec, idx2_vec, val_vec, Nt*ii, prod(obj.resolution)*ii);

        end
        
        function [accelInterpMatrix] = calc_accelInterpmatrix(...
            obj, w000, w001, w010, w011, w100, w101, w110, w111, ii) 
%% Returns the Acceleration Interpolation Matrix
% Calculates and returns the matrix to interpolate the E-Field to get the
% acceleration vector at the particle positions  
            Nt = obj.ParticleArray.Nt;
            
            i000 = sub2ind(obj.resolution, obj.i_x, obj.i_y, obj.i_z);
            i001 = sub2ind(obj.resolution, obj.i_x+1, obj.i_y, obj.i_z);
            i010 = sub2ind(obj.resolution, obj.i_x, obj.i_y+1, obj.i_z);
            i011 = sub2ind(obj.resolution, obj.i_x+1, obj.i_y+1, obj.i_z);
            i100 = sub2ind(obj.resolution, obj.i_x, obj.i_y, obj.i_z+1);
            i101 = sub2ind(obj.resolution, obj.i_x+1, obj.i_y, obj.i_z+1);
            i110 = sub2ind(obj.resolution, obj.i_x, obj.i_y+1, obj.i_z+1);
            i111 = sub2ind(obj.resolution, obj.i_x+1, obj.i_y+1, obj.i_z+1);

            nn = 1:Nt;
%             idx1_vec = [];
%             idx2_vec = [];
%             val_vec = [];
            
           % for ii = 1:obj.ParticleArray.N_particles
                
           idx1 = [nn, nn, nn, nn, nn, nn, nn, nn];
           idx2 = [i000(:,ii); i001(:,ii); i010(:,ii); i011(:,ii); i100(:,ii); i101(:,ii); i110(:,ii); i111(:,ii)];
           val  = [w000(:,ii); w001(:,ii); w010(:,ii); w011(:,ii); w100(:,ii); w101(:,ii); w110(:,ii); w111(:,ii)];
               
               %idx1_vec = [idx1_vec, idx1(:)'];
               %idx2_vec = [idx2_vec, idx2(:)'];
               %val_vec  = [val_vec, val(:)'];
            
           % end
            
            accelInterpMatrix = sparse(idx1, idx2, val, Nt, prod(obj.resolution));

        end
        

        function [D_I_1, D_I_2, D_I_3] = get_I(obj)
%% Calculate Ix, Iy, Iz (Derivative of Trilinear Weights)
% Calculates the derivative of the trilinear weights with respect to x,y,z

            %xv = obj.ParticleArray(1,ii).xv;
            %Nt = obj.ParticleArray(1,ii).Nt;
            
            [D_I_1] = obj.trilinear_weights_I_x;
            [D_I_2] = obj.trilinear_weights_I_y;
            [D_I_3] = obj.trilinear_weights_I_z;

        end 
        
        function obj = set_logvec(obj)
            
            x = obj.ParticleArray.xx;
            y = obj.ParticleArray.yy;
            z = obj.ParticleArray.zz;
            Nt = obj.ParticleArray.Nt;
%             assert(isvector(x))
%             assert(isvector(y))
%             assert(isvector(z))
%             assert(isvector(obj.x_grid))
%             assert(isvector(obj.r_grid))
%             assert(isvector(obj.z_grid))
            assert(isequal(size(x), size(y), size(z)))
            % Coerce x_grid etc. to be row vectors if x, y, z are row vectors...
            if isrow(x) && ~isrow(obj.x_grid)
                row = @(A) reshape(A, 1, []);
                obj.x_grid = row(obj.x_grid);
                obj.r_grid = row(obj.r_grid);
                obj.z_grid = row(obj.z_grid);
            elseif iscolumn(x) && ~iscolumn(obj.x_grid)
                col = @(A) reshape(A, [], 1);
                obj.x_grid = col(obj.x_grid);
                obj.r_grid = col(obj.r_grid);
                obj.z_grid = col(obj.z_grid);
            end


               
            if max(max(x)) >= obj.x_grid(end) || min(min(x)) <= obj.x_grid(1)
                obj.i_x = ones(size(x));
                obj.i_y = ones(size(y));
                obj.i_z = ones(size(z));
                
                    if max(max(y)) >= obj.r_grid(end) || min(min(y)) <= obj.r_grid(1)
                            if max(max(z)) >= obj.z_grid(end) || min(min(z)) <= obj.z_grid(1)
                                  
                                    obj.log_vec = ((x >= obj.x_grid(1)) & (x <= obj.x_grid(end)) & (y >= obj.r_grid(1)) & (y <= obj.r_grid(end)) & (z >= obj.z_grid(1)) & (z <= obj.z_grid(end)));
                                    
                                    obj.i_x(obj.log_vec) = floor( (x(obj.log_vec) - obj.x_grid(1))/obj.dx) +1;
                                    obj.i_y(obj.log_vec) = floor( (y(obj.log_vec) - obj.r_grid(1))/obj.dy) +1;
                                    obj.i_z(obj.log_vec) = floor( (z(obj.log_vec) - obj.z_grid(1))/obj.dz) +1;
                                    
                            else
                                    obj.log_vec = ((x >= obj.x_grid(1)) & (x <= obj.x_grid(end)) & (y >= obj.r_grid(1)) & (y <= obj.r_grid(end)));
                                    
                                    obj.i_x(obj.log_vec) = floor( (x(obj.log_vec) - obj.x_grid(1))/obj.dx) +1;
                                    obj.i_y(obj.log_vec) = floor( (y(obj.log_vec) - obj.r_grid(1))/obj.dy) +1;
                                    obj.i_z(obj.log_vec) = floor( (z(obj.log_vec) - obj.z_grid(1))/obj.dz) +1;
                            end
                            
                    else
                        obj.log_vec = ((x >= obj.x_grid(1)) & (x <= obj.x_grid(end)));

                        obj.i_x(obj.log_vec) = floor( (x(obj.log_vec) - obj.x_grid(1))/obj.dx) +1;
                        obj.i_y(obj.log_vec) = floor( (y(obj.log_vec) - obj.r_grid(1))/obj.dy) +1;
                        obj.i_z(obj.log_vec) = floor( (z(obj.log_vec) - obj.z_grid(1))/obj.dz) +1;
                    end
                        
            else
                
                obj.i_x = floor( (x - obj.x_grid(1))/obj.dx) +1;
                obj.i_y = floor( (y - obj.r_grid(1))/obj.dy) +1;
                obj.i_z = floor( (z - obj.z_grid(1))/obj.dz) +1;
                obj.log_vec = ones(size(obj.i_x));
            end
            
            %log_vec = ((x >= obj.x_grid(1)) & (x <= obj.x_grid(end)) & (y >= obj.y_grid(1)) & (y <= obj.y_grid(end)) & (z >= obj.z_grid(1)) & (z <= obj.z_grid(end)));



            obj.i_x(obj.log_vec) = floor( (x(obj.log_vec) - obj.x_grid(1))/obj.dx) +1;
            obj.i_y(obj.log_vec) = floor( (y(obj.log_vec) - obj.r_grid(1))/obj.dy) +1;
            obj.i_z(obj.log_vec) = floor( (z(obj.log_vec) - obj.z_grid(1))/obj.dz) +1;


            % Handle a special case: when x == obj.x_grid(end) we round its position DOWN
            % instead of UP, to allow all boundary values to be defined as one might
            % expect.  Now x == obj.x_grid(1) is in the first cell AND x == obj.x_grid(end)
            % is in the last cell.

            obj.i_x(x == obj.x_grid(end)) = length(obj.x_grid)-1;
            obj.i_y(y == obj.r_grid(end)) = length(obj.r_grid)-1;
            obj.i_z(z == obj.z_grid(end)) = length(obj.z_grid)-1;

            obj.x_c = obj.x_grid(obj.i_x);
            obj.y_c = obj.r_grid(obj.i_y);
            obj.z_c = obj.z_grid(obj.i_z);
            
        end 
        
        function [D_I_1] = trilinear_weights_I_x(obj)
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

% 
            x = obj.ParticleArray.xx;
            y = obj.ParticleArray.yy;
            z = obj.ParticleArray.zz;


            % Recover weights

            wx = ( 1 ) ./ obj.dx;
            wy = ( y - obj.y_c ) ./ obj.dy;
            wz = ( z - obj.z_c ) ./ obj.dz; 

            w000 = (-wx).*(1-wy).*(1-wz);
            w100 = wx.*(1-wy).*(1-wz);
            w010 = (-wx).*wy.*(1-wz);
            w110 = wx.*wy.*(1-wz);
            w001 = (-wx).*(1-wy).*wz;
            w101 = wx.*(1-wy).*wz;
            w011 = (-wx).*wy.*wz;
            w111 = wx.*wy.*wz;

            w000(~obj.log_vec) = 0;
            w100(~obj.log_vec) = 0;
            w010(~obj.log_vec) = 0;
            w110(~obj.log_vec) = 0;

            w001(~obj.log_vec) = 0;
            w101(~obj.log_vec) = 0;
            w011(~obj.log_vec) = 0;
            w111(~obj.log_vec) = 0;

            assert(isequal(size(obj.i_x), size(x), size(obj.i_y), size(y), size(obj.i_z), size(z), ...
                size(w000), size(w001), size(w010), size(w011), ...
                size(w100), size(w101), size(w110), size(w111)));



           
            i000 = sub2ind(obj.resolution, obj.i_x, obj.i_y, obj.i_z);
            i001 = sub2ind(obj.resolution, obj.i_x, obj.i_y, obj.i_z+1);
            i010 = sub2ind(obj.resolution, obj.i_x, obj.i_y+1, obj.i_z);
            i011 = sub2ind(obj.resolution, obj.i_x, obj.i_y+1, obj.i_z+1);
            i100 = sub2ind(obj.resolution, obj.i_x+1, obj.i_y, obj.i_z);
            i101 = sub2ind(obj.resolution, obj.i_x+1, obj.i_y, obj.i_z+1);
            i110 = sub2ind(obj.resolution, obj.i_x+1, obj.i_y+1, obj.i_z);
            i111 = sub2ind(obj.resolution, obj.i_x+1, obj.i_y+1, obj.i_z+1);

            Nt = obj.ParticleArray.Nt;
            nn = 1:obj.ParticleArray.Nt;
            res = obj.resolution;
            D_I_1 = [];
%             Exx = obj.Ex;
%             Err = obj.Er;
%             Ezz = obj.Ez;

            
            
%             D_I_1_21 = zeros(Nt,1,obj.ParticleArray.N_particles);
%             D_I_1_22 = D_I_1_21;
%             D_I_1_23 = D_I_1_21;
            
            for ii = 1:obj.ParticleArray.N_particles
                
                partial_x_accelInterpMatrix = sparse([nn, nn, nn, nn, nn, nn, nn, nn], ...
                    [i000(:,ii); i001(:,ii); i010(:,ii); i011(:,ii); i100(:,ii); i101(:,ii); i110(:,ii); i111(:,ii)], ...
                    [w000(:,ii); w001(:,ii); w010(:,ii); w011(:,ii); w100(:,ii); w101(:,ii); w110(:,ii); w111(:,ii)], ...
                    Nt, prod(res));
                
               D_I_1 = [D_I_1; partial_x_accelInterpMatrix * obj.Ex(:); partial_x_accelInterpMatrix * obj.Er(:); partial_x_accelInterpMatrix * obj.Ez(:)];
%                D_I_1_21(:,1,ii) = partial_x_accelInterpMatrix * Exx(:);
%                D_I_1_22(:,1,ii) = partial_x_accelInterpMatrix * Err(:);
%                D_I_1_23(:,1,ii) = partial_x_accelInterpMatrix * Ezz(:);

               
            end

            %D_I = cat(2, D_I_1_21, D_I_1_22, D_I_1_23);
            %D_I_1 = D_I(:);
           % D_I_1 = D_I_1_2(:);           

        end     
        

function [D_I_2] = trilinear_weights_I_y(obj)
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
            

            x = obj.ParticleArray.xx;
            y = obj.ParticleArray.yy;
            z = obj.ParticleArray.zz;

            % Recover weights

            wx = ( x - obj.x_c ) ./ obj.dx;
            wy = ( 1 ) ./ obj.dy;
            wz = ( z - obj.z_c ) ./ obj.dz; 

            w000 = (1-wx).*(-wy).*(1-wz);
            w100 = wx.*(-wy).*(1-wz);
            w010 = (1-wx).*wy.*(1-wz);
            w110 = wx.*wy.*(1-wz);
            w001 = (1-wx).*(-wy).*wz;
            w101 = wx.*(-wy).*wz;
            w011 = (1-wx).*wy.*wz;
            w111 = wx.*wy.*wz;

            w000(~obj.log_vec) = 0;
            w100(~obj.log_vec) = 0;
            w010(~obj.log_vec) = 0;
            w110(~obj.log_vec) = 0;

            w001(~obj.log_vec) = 0;
            w101(~obj.log_vec) = 0;
            w011(~obj.log_vec) = 0;
            w111(~obj.log_vec) = 0;

            assert(isequal(size(obj.i_x), size(x), size(obj.i_y), size(y), size(obj.i_z), size(z), ...
                size(w000), size(w001), size(w010), size(w011), ...
                size(w100), size(w101), size(w110), size(w111)));

      
            i000 = sub2ind(obj.resolution, obj.i_x, obj.i_y, obj.i_z);
            i001 = sub2ind(obj.resolution, obj.i_x, obj.i_y, obj.i_z+1);
            i010 = sub2ind(obj.resolution, obj.i_x, obj.i_y+1, obj.i_z);
            i011 = sub2ind(obj.resolution, obj.i_x, obj.i_y+1, obj.i_z+1);
            i100 = sub2ind(obj.resolution, obj.i_x+1, obj.i_y, obj.i_z);
            i101 = sub2ind(obj.resolution, obj.i_x+1, obj.i_y, obj.i_z+1);
            i110 = sub2ind(obj.resolution, obj.i_x+1, obj.i_y+1, obj.i_z);
            i111 = sub2ind(obj.resolution, obj.i_x+1, obj.i_y+1, obj.i_z+1);


            % This Part does not work for Nt = 2
            Nt = obj.ParticleArray.Nt;
            nn = 1:obj.ParticleArray.Nt;
            idx1_vec = [];
            idx2_vec = [];
            val_vec = [];
            
            D_I_2 = [];
            
            for ii = 1:obj.ParticleArray.N_particles
                                
               partial_y_accelInterpMatrix = sparse([nn, nn, nn, nn, nn, nn, nn, nn], ...
                    [i000(:,ii); i001(:,ii); i010(:,ii); i011(:,ii); i100(:,ii); i101(:,ii); i110(:,ii); i111(:,ii)], ...
                    [w000(:,ii); w001(:,ii); w010(:,ii); w011(:,ii); w100(:,ii); w101(:,ii); w110(:,ii); w111(:,ii)], ...
                    obj.ParticleArray.Nt, prod(obj.resolution)); 
                
               D_I_2 = [D_I_2; partial_y_accelInterpMatrix * obj.Ex(:); partial_y_accelInterpMatrix * obj.Er(:); partial_y_accelInterpMatrix * obj.Ez(:)];
      
            end
            

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
        

        function [D_I_3] = trilinear_weights_I_z(obj)
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
            x = obj.ParticleArray.xx;
            y = obj.ParticleArray.yy;
            z = obj.ParticleArray.zz;


            % Recover weights

            wx = ( x - obj.x_c ) ./ obj.dx;
            wy = ( y - obj.y_c ) ./ obj.dy;
            wz = ( 1 ) ./ obj.dz; 

            w000 = (1-wx).*(1-wy).*(-wz);
            w100 = wx.*(1-wy).*(-wz);
            w010 = (1-wx).*wy.*(-wz);
            w110 = wx.*wy.*(-wz);
            w001 = (1-wx).*(1-wy).*wz;
            w101 = wx.*(1-wy).*wz;
            w011 = (1-wx).*wy.*wz;
            w111 = wx.*wy.*wz;

            w000(~obj.log_vec) = 0;
            w100(~obj.log_vec) = 0;
            w010(~obj.log_vec) = 0;
            w110(~obj.log_vec) = 0;

            w001(~obj.log_vec) = 0;
            w101(~obj.log_vec) = 0;
            w011(~obj.log_vec) = 0;
            w111(~obj.log_vec) = 0;

            assert(isequal(size(obj.i_x), size(x), size(obj.i_y), size(y), size(obj.i_z), size(z), ...
                size(w000), size(w001), size(w010), size(w011), ...
                size(w100), size(w101), size(w110), size(w111)));

            i000 = sub2ind(obj.resolution, obj.i_x, obj.i_y, obj.i_z);
            i001 = sub2ind(obj.resolution, obj.i_x, obj.i_y, obj.i_z+1);
            i010 = sub2ind(obj.resolution, obj.i_x, obj.i_y+1, obj.i_z);
            i011 = sub2ind(obj.resolution, obj.i_x, obj.i_y+1, obj.i_z+1);
            i100 = sub2ind(obj.resolution, obj.i_x+1, obj.i_y, obj.i_z);
            i101 = sub2ind(obj.resolution, obj.i_x+1, obj.i_y, obj.i_z+1);
            i110 = sub2ind(obj.resolution, obj.i_x+1, obj.i_y+1, obj.i_z);
            i111 = sub2ind(obj.resolution, obj.i_x+1, obj.i_y+1, obj.i_z+1);


            % This Part does not work for Nt = 2
            Nt = obj.ParticleArray.Nt;
            nn = 1:obj.ParticleArray.Nt;
            
            D_I_3 = [];
            
            for ii = 1:obj.ParticleArray.N_particles               
                
                partial_z_accelInterpMatrix = sparse([nn, nn, nn, nn, nn, nn, nn, nn], ...
                    [i000(:,ii); i001(:,ii); i010(:,ii); i011(:,ii); i100(:,ii); i101(:,ii); i110(:,ii); i111(:,ii)], ...
                    [w000(:,ii); w001(:,ii); w010(:,ii); w011(:,ii); w100(:,ii); w101(:,ii); w110(:,ii); w111(:,ii)], ...
                    obj.ParticleArray.Nt, prod(obj.resolution));
                
                D_I_3 = [D_I_3; partial_z_accelInterpMatrix * obj.Ex(:); partial_z_accelInterpMatrix * obj.Er(:); partial_z_accelInterpMatrix * obj.Ez(:)];
                
                
            end
            
            

        end

        function [systemMatrix, A0, B] = velocityVerletMatrices3D(obj)
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

            [ix_x, ix_y, ix_z, iv_x, iv_y, iv_z] = get_Index3D(obj.ParticleArray.Nt);

            ia_x = ix_x;
            ia_y = ix_y;
            ia_z = ix_z;


            unos = @(n) ones(1,n);

            % ------- Diagonal
    

            ii_diag = [ix_x(1:obj.ParticleArray.Nt), iv_x(1:obj.ParticleArray.Nt)];
            jj_diag = [ix_x(1:obj.ParticleArray.Nt), iv_x(1:obj.ParticleArray.Nt)];
            vv_diag = [-unos(obj.ParticleArray.Nt), -unos(obj.ParticleArray.Nt)];



            ii_A = [ix_x(2:obj.ParticleArray.Nt), ix_x(2:obj.ParticleArray.Nt), iv_x(2:obj.ParticleArray.Nt)];
            jj_A = [ix_x(1:obj.ParticleArray.Nt-1), iv_x(1:obj.ParticleArray.Nt-1), iv_x(1:obj.ParticleArray.Nt-1)];
            vv_A = [unos(obj.ParticleArray.Nt-1), obj.ParticleArray.dt*unos(obj.ParticleArray.Nt-1), unos(obj.ParticleArray.Nt-1)];



            numRows = 6*obj.ParticleArray.Nt;
            numCols = 6*obj.ParticleArray.Nt;

           
            
            idx1_vec = [];
            idx2_vec = [];
            val_vec  = [];
            
            for ii = 1:obj.ParticleArray.N_particles
                
                idx1 = [ii_A, ii_diag, ii_A+obj.ParticleArray.Nt, ii_diag+obj.ParticleArray.Nt, ii_A+2*obj.ParticleArray.Nt, ii_diag+2*obj.ParticleArray.Nt] + numRows*(ii-1);
                idx2 = [jj_A, jj_diag, jj_A+obj.ParticleArray.Nt, jj_diag+obj.ParticleArray.Nt, jj_A+2*obj.ParticleArray.Nt, jj_diag+2*obj.ParticleArray.Nt] + numCols*(ii-1);
                val  = [vv_A, vv_diag, vv_A, vv_diag, vv_A, vv_diag];
                
                idx1_vec = [idx1_vec, idx1(:)'];
                idx2_vec = [idx2_vec, idx2(:)'];
                val_vec  = [val_vec, val(:)'];
                
            end 
            
            systemMatrix = sparse(idx1_vec, idx2_vec, val_vec, ii*numRows, ii*numCols);


            % Right-hand side


            ii_B = [ix_x(2:obj.ParticleArray.Nt), iv_x(2:obj.ParticleArray.Nt), iv_x(2:obj.ParticleArray.Nt)];
            jj_B = [ia_x(1:obj.ParticleArray.Nt-1), ia_x(1:obj.ParticleArray.Nt-1), ia_x(2:obj.ParticleArray.Nt)];
            vv_B = [0.5*obj.ParticleArray.dt*obj.ParticleArray.dt*unos(obj.ParticleArray.Nt-1), 0.5*obj.ParticleArray.dt*unos(obj.ParticleArray.Nt-1), 0.5*obj.ParticleArray.dt*unos(obj.ParticleArray.Nt-1)];

            numRows = 6*obj.ParticleArray.Nt;
            numCols = 3*obj.ParticleArray.Nt;
           
            
            idx1_vec = [];
            idx2_vec = [];
            val_vec  = [];
            
            for ii = 1:obj.ParticleArray.N_particles
                
                idx1 = [ii_B, ii_B+obj.ParticleArray.Nt, ii_B+2*obj.ParticleArray.Nt] + numRows*(ii-1);
                idx2 = [jj_B, jj_B+obj.ParticleArray.Nt, jj_B+2*obj.ParticleArray.Nt] + numCols*(ii-1);
                val  = [vv_B, vv_B, vv_B];
                
                idx1_vec = [idx1_vec, idx1(:)'];
                idx2_vec = [idx2_vec, idx2(:)'];
                val_vec = [val_vec, val(:)'];
                
            end 
            
            B = sparse(idx1_vec, idx2_vec, val_vec, ii*numRows, ii*numCols);

            
            %figure(3); clf
            %imagesc(B)

            % ------- Initial conditions
            % The first timestep needs to take initial conditions:
            % [x(1); v(1)] = A*[x(0); v(0)] + B*[a(0); a(1)]
            % So we need an RHS vector of the right size, 2*Nt x 1.
            % Let's get RHS = A0 * [x(0); v(0)], with A0 being 2*Nt x 2.

            numRows = 6*obj.ParticleArray.Nt;
            numCols = 6;

            % A = [1, dt; 0, 1]
            ii_A0 = [ix_x(1), iv_x(1)];
            jj_A0 = [1, 4];
            vv_A0 = [1, 1];

           
            idx1_vec = [];
            idx2_vec = [];
            val_vec  = [];
            
            for ii = 1:obj.ParticleArray.N_particles
                
                idx1 = [ii_A0, ii_A0+obj.ParticleArray.Nt, ii_A0+2*obj.ParticleArray.Nt] + numRows*(ii-1);
                idx2 = [jj_A0, jj_A0+1, jj_A0+2] + numCols*(ii-1);
                val  = [vv_A0, vv_A0, vv_A0];
                
                idx1_vec = [idx1_vec, idx1(:)'];
                idx2_vec = [idx2_vec, idx2(:)'];
                val_vec = [val_vec, val(:)'];
                
            end 

            A0 = sparse(idx1_vec, idx2_vec, val_vec, ii*numRows, ii*numCols);


        end
        

        
        
        function [obj, S_p] = get_PrimalS2(obj, D_I_1, D_I_2, D_I_3, systemMatrix)
%% Calculation of the Primal Matrix
% This function calculates the matrix of the Primal System (for a
% derivation with respect to the electric field).


               Bp = B_matrix_primal_Array(obj.ParticleArray.t_vec, D_I_1, D_I_2, D_I_3, obj.ParticleArray.N_particles);


               S_p = systemMatrix + Bp;
                
                %end
        end

        
        function [dGdEi_cell] = getdGdE(obj, accelMatrix, ...
             xv_dual)
%% dGdEx, dGdEy, dGdEz Calculation
% This function calculates the derivatives of the objective function
% towards the components of the electric field at each grid point

% 
%             [w000, w001, w010, w011, w100, w101, w110, w111] ...
%                 = trilinear_weights(obj);

           % Nt = obj.ParticleArray.Nt;
            
%             i000 = sub2ind(obj.resolution, obj.i_x, obj.i_y, obj.i_z);
%             i001 = sub2ind(obj.resolution, obj.i_x+1, obj.i_y, obj.i_z);
%             i010 = sub2ind(obj.resolution, obj.i_x, obj.i_y+1, obj.i_z);
%             i011 = sub2ind(obj.resolution, obj.i_x+1, obj.i_y+1, obj.i_z);
%             i100 = sub2ind(obj.resolution, obj.i_x, obj.i_y, obj.i_z+1);
%             i101 = sub2ind(obj.resolution, obj.i_x+1, obj.i_y, obj.i_z+1);
%             i110 = sub2ind(obj.resolution, obj.i_x, obj.i_y+1, obj.i_z+1);
%             i111 = sub2ind(obj.resolution, obj.i_x+1, obj.i_y+1, obj.i_z+1);
% 
%             nn = 1:Nt;
% %             idx1_vec = [];
% %             idx2_vec = [];
% %             val_vec = [];
%             
%            % for ii = 1:obj.ParticleArray.N_particles
%                 
%            idx1 = [nn, nn, nn, nn, nn, nn, nn, nn];
%            idx2 = [i000(:,ii); i001(:,ii); i010(:,ii); i011(:,ii); i100(:,ii); i101(:,ii); i110(:,ii); i111(:,ii)];
%            val  = [w000(:,ii); w001(:,ii); w010(:,ii); w011(:,ii); w100(:,ii); w101(:,ii); w110(:,ii); w111(:,ii)];
%                
%                %idx1_vec = [idx1_vec, idx1(:)'];
%                %idx2_vec = [idx2_vec, idx2(:)'];
%                %val_vec  = [val_vec, val(:)'];
%             
%            % end
%             
%             accelInterpMatrix = sparse(idx1, idx2, val, Nt, prod(obj.resolution));
            
            %end

            dGda = xv_dual' * (-accelMatrix);
            
            [ix_x, ix_y, ix_z, ~] = obj.ParticleArray.get_Index3D();
            
            mult = 1:obj.ParticleArray.N_particles;
            mult = mult * 3* obj.ParticleArray.Nt;
            mult = [0, mult(1:(end-1))];
            idx_x = [];
            idx_y = [];
            idx_z = [];
%             
%             dGda_cell1 = cell(1,obj.ParticleArray.N_particles); 
%             dGda_cell2 = dGda_cell1;
%             dGda_cell3 = dGda_cell2;
%             dGdEi_cell = cell(3,obj.ParticleArray.N_particles); 
            
            for ii = 1:length(mult)
                
%                 dGda_cell1{ii} = dGda(ix_x + mult(ii));
%                 dGda_cell2{ii} = dGda(ix_y + mult(ii));
%                 dGda_cell3{ii} = dGda(ix_z + mult(ii));
                
                idx_x = [idx_x, ix_x + mult(ii)];
                idx_y = [idx_y, ix_y + mult(ii)];
                idx_z = [idx_z, ix_z + mult(ii)];
                
            end
            dGda_cell = cell(1,3);
            dGdEi_cell = cell(1,3);
 
           dGda_cell{1} = dGda(idx_x);
           dGda_cell{2} = dGda(idx_y);
           dGda_cell{3} = dGda(idx_z);
           
           ix = obj.i_x;
           iy = obj.i_y;
           iz = obj.i_z;
           res = obj.resolution;
           xc = obj.x_c;
           yc = obj.y_c;
           zc = obj.z_c;
           x = obj.ParticleArray.xx;
           y = obj.ParticleArray.yy;
           z = obj.ParticleArray.zz;
           ddx = obj.dx;
           ddy = obj.dy;
           ddz = obj.dz;
           llog = obj.log_vec;
           
           
%            i000 = sub2ind(res, ix, iy, iz);
%             i001 = sub2ind(res, ix+1, iy, iz);
%             i010 = sub2ind(res, ix, iy+1, iz);
%             i011 = sub2ind(res, ix+1, iy+1, iz);
%             i100 = sub2ind(res, ix, iy, iz+1);
%             i101 = sub2ind(res, ix+1, iy, iz+1);
%             i110 = sub2ind(res, ix, iy+1, iz+1);
%             i111 = sub2ind(res, ix+1, iy+1, iz+1);

            Nt = obj.ParticleArray.Nt;
            Np = obj.ParticleArray.N_particles;
           
            
                
                


                % Recover weights

                wx = ( x - xc ) ./ ddx;
                wy = ( y - yc ) ./ ddy;
                wz = ( z - zc ) ./ ddz; 

                w000 = (1-wx).*(1-wy).*(1-wz);
                w100 = wx.*(1-wy).*(1-wz);
                w010 = (1-wx).*wy.*(1-wz);
                w110 = wx.*wy.*(1-wz);
                w001 = (1-wx).*(1-wy).*wz;
                w101 = wx.*(1-wy).*wz;
                w011 = (1-wx).*wy.*wz;
                w111 = wx.*wy.*wz;

                w000(~llog) = 0;
                w100(~llog) = 0;
                w010(~llog) = 0;
                w110(~llog) = 0;

                w001(~llog) = 0;
                w101(~llog) = 0;
                w011(~llog) = 0;
                w111(~llog) = 0;
                
            
                i000 = sub2ind(res, ix, iy, iz);
                i001 = sub2ind(res, ix+1, iy, iz);
                i010 = sub2ind(res, ix, iy+1, iz);
                i011 = sub2ind(res, ix+1, iy+1, iz);
                i100 = sub2ind(res, ix, iy, iz+1);
                i101 = sub2ind(res, ix+1, iy, iz+1);
                i110 = sub2ind(res, ix, iy+1, iz+1);
                i111 = sub2ind(res, ix+1, iy+1, iz+1);

                nn = 1:Nt;
                idx1_vec = [];
                idx2_vec = [];
                val_vec = [];

                for ii = 1:Np

                   idx1 = [nn, nn, nn, nn, nn, nn, nn, nn] + Nt*(ii-1);
                   idx2 = [i000(:,ii); i001(:,ii); i010(:,ii); i011(:,ii); i100(:,ii); i101(:,ii); i110(:,ii); i111(:,ii)] + prod(res)*(ii-1);
                   val  = [w000(:,ii); w001(:,ii); w010(:,ii); w011(:,ii); w100(:,ii); w101(:,ii); w110(:,ii); w111(:,ii)];

                   idx1_vec = [idx1_vec, idx1(:)'];
                   idx2_vec = [idx2_vec, idx2(:)'];
                   val_vec  = [val_vec, val(:)'];

                end

                accelInterpMatrix = sparse(idx1_vec, idx2_vec, val_vec, Nt*ii, prod(res)*ii);
                
           for kk = 1:3     
                dGdEi_vec   = dGda_cell{kk} * accelInterpMatrix;
                dGdEi       = sum(reshape(dGdEi_vec, prod(res),...
                    []), 2);
                
                dGdEi_cell{kk}       = reshape(dGdEi, res);
            end
            
            
%             parfor ii = 1:obj.ParticleArray.N_particles
%                 
%                 
% 
%                 nn = 1:Nt;
% %             idx1_vec = [];
% %             idx2_vec = [];
% %             val_vec = [];
%             
%            % for ii = 1:obj.ParticleArray.N_particles
%                 
%             idx1 = [nn, nn, nn, nn, nn, nn, nn, nn];
%             idx2 = [i000(:,ii); i001(:,ii); i010(:,ii); i011(:,ii); i100(:,ii); i101(:,ii); i110(:,ii); i111(:,ii)];
%             val  = [w000(:,ii); w001(:,ii); w010(:,ii); w011(:,ii); w100(:,ii); w101(:,ii); w110(:,ii); w111(:,ii)];
%                
%                %idx1_vec = [idx1_vec, idx1(:)'];
%                %idx2_vec = [idx2_vec, idx2(:)'];
%                %val_vec  = [val_vec, val(:)'];
%             
%            % end
%             
%             accelInterpMatrix = sparse(idx1, idx2, val, Nt, prod(res));
%                 
%                 dGdEx_vec(:,ii) = dGda_cell1{ii} * accelInterpMatrix;
%                 dGdEy_vec(:,ii) = dGda_cell2{ii} * accelInterpMatrix;
%                 dGdEz_vec(:,ii) = dGda_cell3{ii} * accelInterpMatrix;
%                 
%                 
%                 
%             end
%             parfor ii = 1:obj.ParticleArray.N_particles
%                 
%                 [accelInterpMatrix] = obj.calc_accelInterpmatrix(w000, w001, w010, w011, w100, w101, w110, w111, ii);
%                 
%                 %dGdEx_vec{ii} = dGda_cell1{ii} * accelInterpMatrix;
%                 dGdEy_vec(:,ii) = dGda_cell2{ii} * accelInterpMatrix;
%                 %dGdEz_vec{ii} = dGda_cell3{ii} * accelInterpMatrix;
%                 
%                 
%                 
%             end 
%             parfor ii = 1:obj.ParticleArray.N_particles
%                 
%                 [accelInterpMatrix] = obj.calc_accelInterpmatrix(w000, w001, w010, w011, w100, w101, w110, w111, ii);
%                 
%                 %dGdEx_vec{ii} = dGda_cell1{ii} * accelInterpMatrix;
%                 %dGdEy_vec{ii} = dGda_cell2{ii} * accelInterpMatrix;
%                 dGdEz_vec(:,ii) = dGda_cell3{ii} * accelInterpMatrix;
%                 
%                 
%                 
%             end 
            
%             dGdEx = sum(reshape(dGdEx_vec, prod(res),...
%                     []), 2);
%             dGdEy = sum(reshape(dGdEy_vec, prod(res),...
%                 []), 2);
%             dGdEz = sum(reshape(dGdEz_vec, prod(res),...
%                 []), 2);
% 
%             dGdEi_cell{1} = reshape(dGdEx, res);
%             dGdEi_cell{2} = reshape(dGdEy, res);
%             dGdEi_cell{3} = reshape(dGdEz, res);
            
%             for ii = 1:3
%                 
%                 dGdEi_vec   = dGda_cell{ii} * accelInterpMatrix;
%                 dGdEi       = sum(reshape(dGdEi_vec, prod(res),...
%                     []), 2);
%                 
%                 dGdEi_cell{ii}       = reshape(dGdEi, res);
%             end
%             
            
           
            
            
%             dGdEy2 = dGda(idx_y) * accelInterpMatrix;
%             dGdEz2 = dGda(idx_z) * accelInterpMatrix;
%             
%             dGdEx = sum(reshape(dGdEx2, prod(obj.resolution), []), 2);
%             dGdEy = sum(reshape(dGdEy2, prod(obj.resolution), []), 2);
%             dGdEz = sum(reshape(dGdEz2, prod(obj.resolution), []), 2);
%                         
%             dGdEx = reshape(dGdEx, obj.resolution);
%             dGdEy = reshape(dGdEy, obj.resolution);
%             dGdEz = reshape(dGdEz, obj.resolution);


        end
        
        function dGdV_xr = reflect_back(obj,dGdV)

            two_r = size(dGdV,2);
           %dGdV_xr = zeros(size(dGdV,1),round(two_r/2),size(dGdV,3));

            dGdV_xr = dGdV(:,(round(two_r/2):end));
            dGdV_xr(:,2:end) = dGdV_xr(:,2:end) + dGdV(:,(floor(two_r/2):-1:1));

        end 
        
        function plotdFdEx(obj)
            
            figure(1010)
            imagesc(obj.x_grid,obj.r_grid, obj.dFdEx(:,:,1)')
            axis xy image
            colorbar
            xlabel('x [m]')
            ylabel('y [m]')
        end
    end
    
end
