classdef cylindricalVelocity_VerletArray_test
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
        
        D_I_1
        D_I_2
        D_I_3
        Bp
        Bp2
        dGda
        dGdEx
        dGdEy
        dGdEz
        dGdEx2
        dGdEy2
        dGdEz2
        DG
        systemMatrix
        SystemMatrix2
        accelMatrix
        accelMatrix2
        D_I_1_2
        D_I_2_2
        D_I_3_2
        S_p
        S_p2
        Ix
        Iy
        Iz
        
    end
    
    methods        

        function obj = cylindricalVelocity_VerletArray_test(xyz, resolution, ParticleArray, Fobj,varargin)
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
            
            
                obj = obj.calculateF();
                
                %G_sum = 0;

                dGdEx_sum = 0*obj.Ex;
                dGdEr_sum = 0*obj.Er;
                dGdEz_sum = 0*obj.Ez;
                nParticle = size(obj.ParticleArray,2);
            %tt = tic;
           % for ii = 1:length(obj.ParticleArray)
                %ttt = tic;
                %obj.accelFunc = obj.accelerationFunction();
                %obj.ParticleArray.accelerationFunction = obj.accelFunc(obj.Ex, obj.Er, obj.Ez);
                %obj.ParticleArray = obj.ParticleArray.runTrajectory();
                %k(ii) = toc(ttt);
                %tic
                obj = obj.set_logvec;
                
                [w000, w001, w010, w011, w100, w101, w110, w111] ...
                    = trilinear_weights(obj);
                    
                [accelInterpMatrix,accelInterpMatrix2] = obj.get_accelInterpmatrix(w000, w001, w010, w011, w100, w101, w110, w111);
                
                [ix_x, ix_y, ix_z, ~] = obj.ParticleArray.get_Index3D();
                
                [Ix, Iy, Iz, Ix2, Iy2, Iz2, D_I_1, D_I_2, D_I_3] = get_I(obj);
                
                [G_sum, DG] = obj.Fobj(obj.ParticleArray.xv);
                
                DG2 = reshape(DG, [], 1);
                
                % Matricies only depent on the time vector and thus, they
                % are the same for every particle
                [systemMatrix, ~, accelMatrix, systemMatrix2, ~, accelMatrix2] = ...
                    obj.velocityVerletMatrices3D;

                
                
                obj.DG = DG;

                systemMatrix = ...
                    systemMatrix.*...
                    (1/((obj.ParticleArray.ncharges*obj.ParticleArray.e)/(obj.ParticleArray.nmass*obj.ParticleArray.u)));
                
                systemMatrix2 = ...
                    systemMatrix2.*...
                    (1/((obj.ParticleArray.ncharges*obj.ParticleArray.e)/(obj.ParticleArray.nmass*obj.ParticleArray.u)));
                
                [obj, S_p] = obj.get_PrimalS(...
                    Ix, Iy, Iz, systemMatrix);
                
                [obj, S_p2] = obj.get_PrimalS2(...
                    D_I_1, D_I_2, D_I_3, systemMatrix2);
                
                t = 4500; 
                for ii = 1:15
                    
                    isequal(obj.D_I_2{ii},D_I_2(((ii-1)*t+1):(ii*t)))
                    
                end
                
                for ii = 1:obj.ParticleArray.N_particles
                    
                    obj.ParticleArray.xv_dual(:,ii) = S_p{ii}' \ DG(:,ii);
                
                end
                
                obj.systemMatrix = systemMatrix;
                obj.SystemMatrix2 = systemMatrix2;
                obj.accelMatrix = accelMatrix;
                obj.accelMatrix2 = accelMatrix2;
                obj.D_I_1_2 = D_I_1;
                obj.D_I_2_2 = D_I_2;
                obj.D_I_3_2 = D_I_3;
                obj.S_p = S_p;
                obj.S_p2 = S_p2;
                obj.Ix = Ix;
                obj.Iy = Iy;
                obj.Iz = Iz;
                
                xv_dual2 = reshape(S_p2' \ DG2, [], obj.ParticleArray.N_particles);
                xv_dual_vec = S_p2' \ DG2;
                
                a = 9000;
                
                a1 = 1500;
                b1 = 28119240;
                
                a2 = 9000;
                b2 = 4500; 
                
                for ii = 1:15
                    %isequal(obj.Bp{ii}, obj.Bp2((1 + (ii-1)*a):(ii*a),(1 + (ii-1)*a):(ii*a)))
                    %isequal(S_p{ii}, S_p2((1 + (ii-1)*a):(ii*a),(1 + (ii-1)*a):(ii*a)))
                   % isequal(systemMatrix, systemMatrix2((1 + (ii-1)*a):(ii*a),(1 + (ii-1)*a):(ii*a)))
                    isequal(accelInterpMatrix{ii}, accelInterpMatrix2((1 + (ii-1)*a1):(ii*a1),(1 + (ii-1)*b1):(ii*b1)))
                  % isequal(accelMatrix, accelMatrix2((1 + (ii-1)*a2):(ii*a2),(1 + (ii-1)*b2):(ii*b2)))
                  
                end 
                
                dGdEx = {};
                dGdEr = {};
                dGdEz = {};
                dGda = {};
                
                for ii = 1:obj.ParticleArray.N_particles
                    [dGdEx{ii}, dGdEr{ii}, dGdEz{ii}, dGda{ii}] = obj.getdGdE(ii, accelMatrix, ix_x, ix_y, ix_z, accelInterpMatrix);
                    dGdEx_sum = dGdEx_sum + dGdEx{ii};
                    dGdEr_sum = dGdEr_sum + dGdEr{ii};
                    dGdEz_sum = dGdEz_sum + dGdEz{ii};
                end
                
                obj = obj.getdGdE2(accelMatrix2, accelInterpMatrix2, xv_dual_vec);
                
              %  obj.dGdEx = obj.dGdEx *1/nParticle*2*G;
               % obj.dGdEy = obj.dGdEy *1/nParticle*2*G;
               % obj.dGdEz = obj.dGdEz *1/nParticle*2*G;
                for ii = 1:15
                  
                   % isequal(dGda{ii}, obj.dGda((1 + (ii-1)*b2):b2*ii))
%                     isequal(dGdEx{ii}(:)', obj.dGdEx2((1 + (ii-1)*b1):b1*ii))
%                     isequal(dGdEr{ii}(:)', obj.dGdEy2((1 + (ii-1)*b1):b1*ii))
%                     isequal(dGdEz{ii}(:)', obj.dGdEz2((1 + (ii-1)*b1):b1*ii))
                    
                    isequal(dGdEx_sum, obj.dGdEx)
                    isequal(dGdEr_sum, obj.dGdEy)
                    isequal(dGdEz_sum, obj.dGdEz)

                end 
                
              %  dGdEx_sum2 = ;
              %  dGdEr_sum2 = ;
              %  dGdEz_sum2 = ;
              
            mult = 1:obj.ParticleArray.N_particles;
            mult = mult * 3* obj.ParticleArray.Nt;
            mult = [0, mult(1:(end-1))];
            idx_x = [];
            idx_y = [];
            idx_z = [];
            
            for ii = 1:length(mult)
                
                idx_x = [idx_x, ix_x + mult(ii)];
                idx_y = [idx_y, ix_y + mult(ii)];
                idx_z = [idx_z, ix_z + mult(ii)];
                
            end
                
                %G_sum = G_sum + G.^2;
        
               % dGdEx_sum = dGdEx_sum + 1/nParticle*2*dGdEx*G;
              %  dGdEr_sum = dGdEr_sum + 1/nParticle*2*dGdEr*G;
               % dGdEz_sum = dGdEz_sum + 1/nParticle*2*dGdEz*G;
                
           % end
           
           dGdEx_sum = obj.dGdEx;
           dGdEr_sum = obj.dGdEy;
           dGdEz_sum = obj.dGdEz;
            
            obj.Fval = G_sum;
            %dGdEx_sum = 0.5*dGdEx_sum ./ obj.Fval;
            %dGdEr_sum = 0.5*dGdEr_sum ./ obj.Fval;
            %dGdEz_sum = 0.5*dGdEz_sum ./ obj.Fval;
            
            
            dGdV = 0*dGdEx_sum;
            
            dGdV(3:end,:,:)   = dGdV(3:end,:,:)     +...
                (-0.5/obj.dx)*dGdEx_sum(2:end-1,:,:);
            dGdV(1:end-2,:,:) = dGdV(1:end-2,:,:)   +...
                (0.5/obj.dx)*dGdEx_sum(2:end-1,:,:);

            dGdV(:,3:end,:)   = dGdV(:,3:end,:)     +...
                (-0.5/obj.dy)*dGdEr_sum(:,2:end-1,:);
            dGdV(:,1:end-2,:) = dGdV(:,1:end-2,:)   +...
                (0.5/obj.dy)*dGdEr_sum(:,2:end-1,:);

            dGdV(:,:,3:end)   = dGdV(:,:,3:end)     +...
                (-0.5/obj.dz)*dGdEz_sum(:,:,2:end-1);
            dGdV(:,:,1:end-2) = dGdV(:,:,1:end-2)   +...
                (0.5/obj.dz)*dGdEz_sum(:,:,2:end-1);
            
           
            dFdEx_full = sum(dGdEx_sum,3);
            dFdEr_full = sum(dGdEr_sum,3);
            dFdEz_full = sum(dGdEz_sum,3);
            
            dFdV_full = sum(dGdV,3);
            
            obj.dFdEx = obj.reflect_back(dFdEx_full);
            obj.dFdEr = obj.reflect_back(dFdEr_full);
            obj.dFdEz = obj.reflect_back(dFdEz_full);
           
            obj.dFdV = obj.reflect_back(dFdV_full);
          %  disp('forward VV time')
           % sum(k)
          %  disp('dual vv time')
           % dual_time = toc(tt) - sum(k);
          %  disp(dual_time)
            
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
        
        
        
        
        function accelFunc = accelerationFunction2(obj)
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
            interpolant(Ex, xyz(1), xyz(2), xyz(3)); ...
            interpolant(Ey, xyz(1), xyz(2), xyz(3)); ...
            interpolant(Ez, xyz(1), xyz(2), xyz(3)); ...
            ];

        end
        
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

            %             assert(isvector(x))
            %             assert(isvector(y))
            %             assert(isvector(z))
            %             assert(isvector(x_grid))
            %             assert(isvector(y_grid))
            %             assert(isvector(z_grid))
            %             assert(isequal(size(x), size(y), size(z)))
            x = obj.ParticleArray.xx;
            y = obj.ParticleArray.yy;
            z = obj.ParticleArray.zz;
%                         % Coerce x_grid etc. to be row vectors if x, y, z are row vectors...
%             if isrow(x) && ~isrow(obj.x_grid)
%                 row = @(A) reshape(A, 1, []);
%                 obj.x_grid = row(obj.x_grid);
%                 obj.r_grid = row(obj.r_grid);
%                 obj.z_grid = row(obj.z_grid);
%             elseif iscolumn(x) && ~iscolumn(obj.x_grid)
%                 col = @(A) reshape(A, [], 1);
%                 obj.x_grid = col(obj.x_grid);
%                 obj.r_grid = col(obj.r_grid);
%                 obj.z_grid = col(obj.z_grid);
%             end
% 
% 
%             %             d_x = x_grid(2)-x_grid(1);
%             %             d_y = y_grid(2)-y_grid(1);
%             %             d_z = z_grid(2)-z_grid(1);
%             
% 
% 
% 
%             if max(max(x)) >= obj.x_grid(end) || min(min(x)) <= obj.x_grid(1)
%                 i_x = ones(size(x));
%                 i_y = ones(size(y));
%                 i_z = ones(size(z));
%                     if max(max(y)) >= obj.r_grid(end) || min(min(y)) <= obj.r_grid(1)
%                             if max(max(z)) >= obj.z_grid(end) || min(min(z)) <= obj.z_grid(1)
%                                   
%                                     log_vec = ((x >= obj.x_grid(1)) & (x <= obj.x_grid(end)) & (y >= obj.r_grid(1)) & (y <= obj.r_grid(end)) & (z >= obj.z_grid(1)) & (z <= obj.z_grid(end)));
%                                     
%                                     i_x(log_vec) = floor( (x(log_vec) - obj.x_grid(1))/obj.dx) +1;
%                                     i_y(log_vec) = floor( (y(log_vec) - obj.r_grid(1))/obj.dy) +1;
%                                     i_z(log_vec) = floor( (z(log_vec) - obj.z_grid(1))/obj.dz) +1;
%                                     
%                             else
%                                     log_vec = ((x >= obj.x_grid(1)) & (x <= obj.x_grid(end)) & (y >= obj.r_grid(1)) & (y <= obj.r_grid(end)));
%                                     
%                                     i_x(log_vec) = floor( (x(log_vec) - obj.x_grid(1))/obj.dx) +1;
%                                     i_y(log_vec) = floor( (y(log_vec) - obj.r_grid(1))/obj.dy) +1;
%                                     i_z(log_vec) = floor( (z(log_vec) - obj.z_grid(1))/obj.dz) +1;
%                             end
%                             
%                     else
%                         log_vec = ((x >= obj.x_grid(1)) & (x <= obj.x_grid(end)));
% 
%                         i_x(log_vec) = floor( (x(log_vec) - obj.x_grid(1))/obj.dx) +1;
%                         i_y(log_vec) = floor( (y(log_vec) - obj.r_grid(1))/obj.dy) +1;
%                         i_z(log_vec) = floor( (z(log_vec) - obj.z_grid(1))/obj.dz) +1;
%                     end
%                         
%             else
%                 
%                 i_x = floor( (x - obj.x_grid(1))/obj.dx) +1;
%                 i_y = floor( (y - obj.r_grid(1))/obj.dy) +1;
%                 i_z = floor( (z - obj.z_grid(1))/obj.dz) +1;
%                 log_vec = ones(size(i_x));
%             end
% 
%             % Handle a special case: when x == x_grid(end) we round its position DOWN
%             % instead of UP, to allow all boundary values to be defined as one might
%             % expect.  Now x == x_grid(1) is in the first cell AND x == x_grid(end)
%             % is in the last cell.
% 
%             i_x(x == obj.x_grid(end)) = length(obj.x_grid)-1;
%             i_y(y == obj.r_grid(end)) = length(obj.r_grid)-1;
%             i_z(z == obj.z_grid(end)) = length(obj.z_grid)-1;
% 
%             x_c = obj.x_grid(i_x);
%             y_c = obj.r_grid(i_y);
%             z_c = obj.z_grid(i_z);

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
        
        function [accelInterpMatrix, accelInterpMatrix2] = get_accelInterpmatrix(...
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
                accelInterpMatrix{ii} = sparse([nn, nn, nn, nn, nn, nn, nn, nn], ...
                    [i000(:,ii); i001(:,ii); i010(:,ii); i011(:,ii); i100(:,ii); i101(:,ii); i110(:,ii); i111(:,ii)], ...
                    [w000(:,ii); w001(:,ii); w010(:,ii); w011(:,ii); w100(:,ii); w101(:,ii); w110(:,ii); w111(:,ii)], ...
                    Nt, prod(obj.resolution));
               idx1 = [nn, nn, nn, nn, nn, nn, nn, nn] + Nt*(ii-1);
               idx2 = [i000(:,ii); i001(:,ii); i010(:,ii); i011(:,ii); i100(:,ii); i101(:,ii); i110(:,ii); i111(:,ii)] + prod(obj.resolution)*(ii-1);
               val  = [w000(:,ii); w001(:,ii); w010(:,ii); w011(:,ii); w100(:,ii); w101(:,ii); w110(:,ii); w111(:,ii)];
               
               idx1_vec = [idx1_vec, idx1(:)'];
               idx2_vec = [idx2_vec, idx2(:)'];
               val_vec  = [val_vec, val(:)'];
            
            end
            
            accelInterpMatrix2 = sparse(idx1_vec, idx2_vec, val_vec, Nt*ii, prod(obj.resolution)*ii);

        end
        

        function [Ix, Iy, Iz, Ix2, Iy2, Iz2, D_I_1, D_I_2, D_I_3] = get_I(obj)
%% Calculate Ix, Iy, Iz (Derivative of Trilinear Weights)
% Calculates the derivative of the trilinear weights with respect to x,y,z

            %xv = obj.ParticleArray(1,ii).xv;
            %Nt = obj.ParticleArray(1,ii).Nt;
            
            [Ix, Ix2, D_I_1] = obj.trilinear_weights_I_x;
            [Iy, Iy2, D_I_2] = obj.trilinear_weights_I_y;
            [Iz, Iz2, D_I_3] = obj.trilinear_weights_I_z;

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
        
        function [partial_x_accelInterpMatrix, partial_x_accelInterpMatrix2, D_I_1, w000, w001, w010, w011, w100, w101, w110, w111] = trilinear_weights_I_x(obj)
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
%             Nt = obj.ParticleArray.Nt;
% %             assert(isvector(x))
% %             assert(isvector(y))
% %             assert(isvector(z))
% %             assert(isvector(obj.x_grid))
% %             assert(isvector(obj.r_grid))
% %             assert(isvector(obj.z_grid))
%             assert(isequal(size(x), size(y), size(z)))
%             % Coerce x_grid etc. to be row vectors if x, y, z are row vectors...
%             if isrow(x) && ~isrow(obj.x_grid)
%                 row = @(A) reshape(A, 1, []);
%                 obj.x_grid = row(obj.x_grid);
%                 obj.r_grid = row(obj.r_grid);
%                 obj.z_grid = row(obj.z_grid);
%             elseif iscolumn(x) && ~iscolumn(obj.x_grid)
%                 col = @(A) reshape(A, [], 1);
%                 obj.x_grid = col(obj.x_grid);
%                 obj.r_grid = col(obj.r_grid);
%                 obj.z_grid = col(obj.z_grid);
%             end
% 
% 
%                
%             if max(max(x)) >= obj.x_grid(end) || min(min(x)) <= obj.x_grid(1)
%                 i_x = ones(size(x));
%                 i_y = ones(size(y));
%                 i_z = ones(size(z));
%                 
%                     if max(max(y)) >= obj.r_grid(end) || min(min(y)) <= obj.r_grid(1)
%                             if max(max(z)) >= obj.z_grid(end) || min(min(z)) <= obj.z_grid(1)
%                                   
%                                     log_vec = ((x >= obj.x_grid(1)) & (x <= obj.x_grid(end)) & (y >= obj.r_grid(1)) & (y <= obj.r_grid(end)) & (z >= obj.z_grid(1)) & (z <= obj.z_grid(end)));
%                                     
%                                     i_x(log_vec) = floor( (x(log_vec) - obj.x_grid(1))/obj.dx) +1;
%                                     i_y(log_vec) = floor( (y(log_vec) - obj.r_grid(1))/obj.dy) +1;
%                                     i_z(log_vec) = floor( (z(log_vec) - obj.z_grid(1))/obj.dz) +1;
%                                     
%                             else
%                                     log_vec = ((x >= obj.x_grid(1)) & (x <= obj.x_grid(end)) & (y >= obj.r_grid(1)) & (y <= obj.r_grid(end)));
%                                     
%                                     i_x(log_vec) = floor( (x(log_vec) - obj.x_grid(1))/obj.dx) +1;
%                                     i_y(log_vec) = floor( (y(log_vec) - obj.r_grid(1))/obj.dy) +1;
%                                     i_z(log_vec) = floor( (z(log_vec) - obj.z_grid(1))/obj.dz) +1;
%                             end
%                             
%                     else
%                         log_vec = ((x >= obj.x_grid(1)) & (x <= obj.x_grid(end)));
% 
%                         i_x(log_vec) = floor( (x(log_vec) - obj.x_grid(1))/obj.dx) +1;
%                         i_y(log_vec) = floor( (y(log_vec) - obj.r_grid(1))/obj.dy) +1;
%                         i_z(log_vec) = floor( (z(log_vec) - obj.z_grid(1))/obj.dz) +1;
%                     end
%                         
%             else
%                 
%                 i_x = floor( (x - obj.x_grid(1))/obj.dx) +1;
%                 i_y = floor( (y - obj.r_grid(1))/obj.dy) +1;
%                 i_z = floor( (z - obj.z_grid(1))/obj.dz) +1;
%                 log_vec = ones(size(i_x));
%             end
%             
%             %log_vec = ((x >= obj.x_grid(1)) & (x <= obj.x_grid(end)) & (y >= obj.y_grid(1)) & (y <= obj.y_grid(end)) & (z >= obj.z_grid(1)) & (z <= obj.z_grid(end)));
% 
% 
% 
%             i_x(log_vec) = floor( (x(log_vec) - obj.x_grid(1))/obj.dx) +1;
%             i_y(log_vec) = floor( (y(log_vec) - obj.r_grid(1))/obj.dy) +1;
%             i_z(log_vec) = floor( (z(log_vec) - obj.z_grid(1))/obj.dz) +1;
% 
% 
%             % Handle a special case: when x == obj.x_grid(end) we round its position DOWN
%             % instead of UP, to allow all boundary values to be defined as one might
%             % expect.  Now x == obj.x_grid(1) is in the first cell AND x == obj.x_grid(end)
%             % is in the last cell.
% 
%             i_x(x == obj.x_grid(end)) = length(obj.x_grid)-1;
%             i_y(y == obj.r_grid(end)) = length(obj.r_grid)-1;
%             i_z(z == obj.z_grid(end)) = length(obj.z_grid)-1;
% 
%             x_c = obj.x_grid(i_x);
%             y_c = obj.r_grid(i_y);
%             z_c = obj.z_grid(i_z);

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
            idx1_vec = [];
            idx2_vec = [];
            val_vec = [];
            
            D_I_1 = [];
            
            for ii = 1:obj.ParticleArray.N_particles
                partial_x_accelInterpMatrix{ii} = sparse([nn, nn, nn, nn, nn, nn, nn, nn], ...
                    [i000(:,ii); i001(:,ii); i010(:,ii); i011(:,ii); i100(:,ii); i101(:,ii); i110(:,ii); i111(:,ii)], ...
                    [w000(:,ii); w001(:,ii); w010(:,ii); w011(:,ii); w100(:,ii); w101(:,ii); w110(:,ii); w111(:,ii)], ...
                    obj.ParticleArray.Nt, prod(obj.resolution));
                
                
                partial_x_accelInterpMatrix3 = sparse([nn, nn, nn, nn, nn, nn, nn, nn], ...
                    [i000(:,ii); i001(:,ii); i010(:,ii); i011(:,ii); i100(:,ii); i101(:,ii); i110(:,ii); i111(:,ii)], ...
                    [w000(:,ii); w001(:,ii); w010(:,ii); w011(:,ii); w100(:,ii); w101(:,ii); w110(:,ii); w111(:,ii)], ...
                    obj.ParticleArray.Nt, prod(obj.resolution));
                
               D_I_1 = [D_I_1; partial_x_accelInterpMatrix3 * obj.Ex(:); partial_x_accelInterpMatrix3 * obj.Er(:); partial_x_accelInterpMatrix3 * obj.Ez(:)];

                
               idx1 = [nn, nn, nn, nn, nn, nn, nn, nn] + Nt*(ii-1);
               idx2 = [i000(:,ii); i001(:,ii); i010(:,ii); i011(:,ii); i100(:,ii); i101(:,ii); i110(:,ii); i111(:,ii)] + prod(obj.resolution)*(ii-1);
               val  = [w000(:,ii); w001(:,ii); w010(:,ii); w011(:,ii); w100(:,ii); w101(:,ii); w110(:,ii); w111(:,ii)];
               
               idx1_vec = [idx1_vec, idx1(:)];
               idx2_vec = [idx2_vec, idx2(:)];
               val_vec  = [val_vec, val(:)];
               
            end
            
            partial_x_accelInterpMatrix2 = sparse(idx1_vec, idx2_vec, val_vec, Nt*ii, prod(obj.resolution)*ii);

        end     
        

function [partial_y_accelInterpMatrix, partial_y_accelInterpMatrix2, D_I_2, w000, w001, w010, w011, w100, w101, w110, w111] = trilinear_weights_I_y(obj)
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
%             Nt = obj.ParticleArray(ii).Nt;
%             
% %             assert(isvector(x))
% %             assert(isvector(y))
% %             assert(isvector(z))
% %             assert(isvector(obj.x_grid))
% %             assert(isvector(obj.r_grid))
% %             assert(isvector(obj.z_grid))
% %             assert(isequal(size(x), size(y), size(z)))
% 
%             % Coerce obj.x_grid etc. to be row vectors if x, y, z are row vectors...
%             if isrow(x) && ~isrow(obj.x_grid)
%                 row = @(A) reshape(A, 1, []);
%                 obj.x_grid = row(obj.x_grid);
%                 obj.r_grid = row(obj.r_grid);
%                 obj.z_grid = row(obj.z_grid);
%             elseif iscolumn(x) && ~iscolumn(obj.x_grid)
%                 col = @(A) reshape(A, [], 1);
%                 obj.x_grid = col(obj.x_grid);
%                 obj.r_grid = col(obj.r_grid);
%                 obj.z_grid = col(obj.z_grid);
%             end
% 
% 
%             % TODO: what does "log" mean in "log_vec"? - It was supposed to mean
%             % "Logic" just because it is a binary vector
%              if max(max(x)) >= obj.x_grid(end) || min(min(x)) <= obj.x_grid(1)
%                 i_x = ones(size(x));
%                 i_y = ones(size(y));
%                 i_z = ones(size(z));
%                     if max(max(y)) >= obj.r_grid(end) || min(min(y)) <= obj.r_grid(1)
%                             if max(max(z)) >= obj.z_grid(end) || min(min(z)) <= obj.z_grid(1)
%                                   
%                                     log_vec = ((x >= obj.x_grid(1)) & (x <= obj.x_grid(end)) & (y >= obj.r_grid(1)) & (y <= obj.r_grid(end)) & (z >= obj.z_grid(1)) & (z <= obj.z_grid(end)));
%                                     
%                                     i_x(log_vec) = floor( (x(log_vec) - obj.x_grid(1))/obj.dx) +1;
%                                     i_y(log_vec) = floor( (y(log_vec) - obj.r_grid(1))/obj.dy) +1;
%                                     i_z(log_vec) = floor( (z(log_vec) - obj.z_grid(1))/obj.dz) +1;
%                                     
%                             else
%                                     log_vec = ((x >= obj.x_grid(1)) & (x <= obj.x_grid(end)) & (y >= obj.r_grid(1)) & (y <= obj.r_grid(end)));
%                                     
%                                     i_x(log_vec) = floor( (x(log_vec) - obj.x_grid(1))/obj.dx) +1;
%                                     i_y(log_vec) = floor( (y(log_vec) - obj.r_grid(1))/obj.dy) +1;
%                                     i_z(log_vec) = floor( (z(log_vec) - obj.z_grid(1))/obj.dz) +1;
%                             end
%                             
%                     else
%                         log_vec = ((x >= obj.x_grid(1)) & (x <= obj.x_grid(end)));
% 
%                         i_x(log_vec) = floor( (x(log_vec) - obj.x_grid(1))/obj.dx) +1;
%                         i_y(log_vec) = floor( (y(log_vec) - obj.r_grid(1))/obj.dy) +1;
%                         i_z(log_vec) = floor( (z(log_vec) - obj.z_grid(1))/obj.dz) +1;
%                     end
%                         
%             else
%                 
%                 i_x = floor( (x - obj.x_grid(1))/obj.dx) +1;
%                 i_y = floor( (y - obj.r_grid(1))/obj.dy) +1;
%                 i_z = floor( (z - obj.z_grid(1))/obj.dz) +1;
%                 log_vec = ones(size(i_x));
% 
%             end
% 
%      
% 
%             % Handle a special case: when x == obj.x_grid(end) we round its position DOWN
%             % instead of UP, to allow all boundary values to be defined as one might
%             % expect.  Now x == obj.x_grid(1) is in the first cell AND x == obj.x_grid(end)
%             % is in the last cell.
% 
%             i_x(x == obj.x_grid(end)) = length(obj.x_grid)-1;
%             i_y(y == obj.r_grid(end)) = length(obj.r_grid)-1;
%             i_z(z == obj.z_grid(end)) = length(obj.z_grid)-1;
% 
%             x_c = obj.x_grid(i_x);
%             y_c = obj.r_grid(i_y);
%             z_c = obj.z_grid(i_z);

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
                partial_y_accelInterpMatrix{ii} = sparse([nn, nn, nn, nn, nn, nn, nn, nn], ...
                    [i000(:,ii); i001(:,ii); i010(:,ii); i011(:,ii); i100(:,ii); i101(:,ii); i110(:,ii); i111(:,ii)], ...
                    [w000(:,ii); w001(:,ii); w010(:,ii); w011(:,ii); w100(:,ii); w101(:,ii); w110(:,ii); w111(:,ii)], ...
                    obj.ParticleArray.Nt, prod(obj.resolution));  
                
               partial_y_accelInterpMatrix3 = sparse([nn, nn, nn, nn, nn, nn, nn, nn], ...
                    [i000(:,ii); i001(:,ii); i010(:,ii); i011(:,ii); i100(:,ii); i101(:,ii); i110(:,ii); i111(:,ii)], ...
                    [w000(:,ii); w001(:,ii); w010(:,ii); w011(:,ii); w100(:,ii); w101(:,ii); w110(:,ii); w111(:,ii)], ...
                    obj.ParticleArray.Nt, prod(obj.resolution)); 
                
               D_I_2 = [D_I_2; partial_y_accelInterpMatrix3 * obj.Ex(:); partial_y_accelInterpMatrix3 * obj.Er(:); partial_y_accelInterpMatrix3 * obj.Ez(:)];

               idx1 = [nn, nn, nn, nn, nn, nn, nn, nn] + Nt*(ii-1);
               idx2 = [i000(:,ii); i001(:,ii); i010(:,ii); i011(:,ii); i100(:,ii); i101(:,ii); i110(:,ii); i111(:,ii)] + prod(obj.resolution)*(ii-1);
               val  = [w000(:,ii); w001(:,ii); w010(:,ii); w011(:,ii); w100(:,ii); w101(:,ii); w110(:,ii); w111(:,ii)];
               
               idx1_vec = [idx1_vec, idx1(:)];
               idx2_vec = [idx2_vec, idx2(:)];
               val_vec  = [val_vec, val(:)];
            end
            
            partial_y_accelInterpMatrix2 = sparse(idx1_vec, idx2_vec, val_vec, Nt*ii, prod(obj.resolution)*ii);

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
        

        function [partial_z_accelInterpMatrix, partial_z_accelInterpMatrix2, D_I_3, w000, w001, w010, w011, w100, w101, w110, w111] = trilinear_weights_I_z(obj)
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
%             Nt = obj.ParticleArray(ii).Nt;
%             
%             assert(isvector(x))
%             assert(isvector(y))
%             assert(isvector(z))
%             assert(isvector(obj.x_grid))
%             assert(isvector(obj.r_grid))
%             assert(isvector(obj.z_grid))
%             assert(isequal(size(x), size(y), size(z)))
% 
%             % Coerce obj.x_grid etc. to be row vectors if x, y, z are row vectors...
%             if isrow(x) && ~isrow(obj.x_grid)
%                 row = @(A) reshape(A, 1, []);
%                 obj.x_grid = row(obj.x_grid);
%                 obj.r_grid = row(obj.r_grid);
%                 obj.z_grid = row(obj.z_grid);
%             elseif iscolumn(x) && ~iscolumn(obj.x_grid)
%                 col = @(A) reshape(A, [], 1);
%                 obj.x_grid = col(obj.x_grid);
%                 obj.r_grid = col(obj.r_grid);
%                 obj.z_grid = col(obj.z_grid);
%             end
% 
% 
% 
%            if max(x) >= obj.x_grid(end) || min(x) <= obj.x_grid(1)
%                 i_x = ones(size(x));
%                 i_y = ones(size(y));
%                 i_z = ones(size(z));
%                     if max(y) >= obj.r_grid(end) || min(y) <= obj.r_grid(1)
%                             if max(z) >= obj.z_grid(end) || min(z) <= obj.z_grid(1)
%                                   
%                                     log_vec = ((x >= obj.x_grid(1)) & (x <= obj.x_grid(end)) & (y >= obj.r_grid(1)) & (y <= obj.r_grid(end)) & (z >= obj.z_grid(1)) & (z <= obj.z_grid(end)));
%                                     
%                                     i_x(log_vec) = floor( (x(log_vec) - obj.x_grid(1))/obj.dx) +1;
%                                     i_y(log_vec) = floor( (y(log_vec) - obj.r_grid(1))/obj.dy) +1;
%                                     i_z(log_vec) = floor( (z(log_vec) - obj.z_grid(1))/obj.dz) +1;
%                                     
%                             else
%                                     log_vec = ((x >= obj.x_grid(1)) & (x <= obj.x_grid(end)) & (y >= obj.r_grid(1)) & (y <= obj.r_grid(end)));
%                                     
%                                     i_x(log_vec) = floor( (x(log_vec) - obj.x_grid(1))/obj.dx) +1;
%                                     i_y(log_vec) = floor( (y(log_vec) - obj.r_grid(1))/obj.dy) +1;
%                                     i_z(log_vec) = floor( (z(log_vec) - obj.z_grid(1))/obj.dz) +1;
%                             end
%                             
%                     else
%                         log_vec = ((x >= obj.x_grid(1)) & (x <= obj.x_grid(end)));
% 
%                         i_x(log_vec) = floor( (x(log_vec) - obj.x_grid(1))/obj.dx) +1;
%                         i_y(log_vec) = floor( (y(log_vec) - obj.r_grid(1))/obj.dy) +1;
%                         i_z(log_vec) = floor( (z(log_vec) - obj.z_grid(1))/obj.dz) +1;
%                     end
%                         
%             else
%                 
%                 i_x = floor( (x - obj.x_grid(1))/obj.dx) +1;
%                 i_y = floor( (y - obj.r_grid(1))/obj.dy) +1;
%                 i_z = floor( (z - obj.z_grid(1))/obj.dz) +1;
%                 log_vec = ones(size(i_x));
%             end
% 
%             % Handle a special case: when x == obj.x_grid(end) we round its position DOWN
%             % instead of UP, to allow all boundary values to be defined as one might
%             % expect.  Now x == obj.x_grid(1) is in the first cell AND x == obj.x_grid(end)
%             % is in the last cell.
% 
%             i_x(x == obj.x_grid(end)) = length(obj.x_grid)-1;
%             i_y(y == obj.r_grid(end)) = length(obj.r_grid)-1;
%             i_z(z == obj.z_grid(end)) = length(obj.z_grid)-1;
% 
%             x_c = obj.x_grid(i_x);
%             y_c = obj.r_grid(i_y);
%             z_c = obj.z_grid(i_z);

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
            idx1_vec = [];
            idx2_vec = [];
            val_vec = [];
            
            D_I_3 = [];
            
            for ii = 1:obj.ParticleArray.N_particles
                partial_z_accelInterpMatrix{ii} = sparse([nn, nn, nn, nn, nn, nn, nn, nn], ...
                    [i000(:,ii); i001(:,ii); i010(:,ii); i011(:,ii); i100(:,ii); i101(:,ii); i110(:,ii); i111(:,ii)], ...
                    [w000(:,ii); w001(:,ii); w010(:,ii); w011(:,ii); w100(:,ii); w101(:,ii); w110(:,ii); w111(:,ii)], ...
                    obj.ParticleArray.Nt, prod(obj.resolution));
                
                
                partial_z_accelInterpMatrix3 = sparse([nn, nn, nn, nn, nn, nn, nn, nn], ...
                    [i000(:,ii); i001(:,ii); i010(:,ii); i011(:,ii); i100(:,ii); i101(:,ii); i110(:,ii); i111(:,ii)], ...
                    [w000(:,ii); w001(:,ii); w010(:,ii); w011(:,ii); w100(:,ii); w101(:,ii); w110(:,ii); w111(:,ii)], ...
                    obj.ParticleArray.Nt, prod(obj.resolution));
                
                D_I_3 = [D_I_3; partial_z_accelInterpMatrix3 * obj.Ex(:); partial_z_accelInterpMatrix3 * obj.Er(:); partial_z_accelInterpMatrix3 * obj.Ez(:)];
                
                idx1 = [nn, nn, nn, nn, nn, nn, nn, nn] + Nt*(ii-1);
                idx2 = [i000(:,ii); i001(:,ii); i010(:,ii); i011(:,ii); i100(:,ii); i101(:,ii); i110(:,ii); i111(:,ii)] + prod(obj.resolution)*(ii-1);
                val  = [w000(:,ii); w001(:,ii); w010(:,ii); w011(:,ii); w100(:,ii); w101(:,ii); w110(:,ii); w111(:,ii)];
               
                idx1_vec = [idx1_vec, idx1(:)];
                idx2_vec = [idx2_vec, idx2(:)];
                val_vec  = [val_vec, val(:)];
            end
            
            partial_z_accelInterpMatrix2 = sparse(idx1_vec, idx2_vec, val_vec, Nt*ii, prod(obj.resolution)*ii);

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

        function [systemMatrix, A0, B, systemMatrix2, A02, B2] = velocityVerletMatrices3D(obj)
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

            systemMatrix = sparse(...
                [ii_A, ii_diag, ii_A+obj.ParticleArray.Nt, ii_diag+obj.ParticleArray.Nt, ii_A+2*obj.ParticleArray.Nt, ii_diag+2*obj.ParticleArray.Nt], ...
                [jj_A, jj_diag, jj_A+obj.ParticleArray.Nt, jj_diag+obj.ParticleArray.Nt, jj_A+2*obj.ParticleArray.Nt, jj_diag+2*obj.ParticleArray.Nt], ...
                [vv_A, vv_diag, vv_A, vv_diag, vv_A, vv_diag],...
                numRows, numCols);
            
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
            
            systemMatrix2 = sparse(idx1_vec, idx2_vec, val_vec, ii*numRows, ii*numCols);


            % Right-hand side


            ii_B = [ix_x(2:obj.ParticleArray.Nt), iv_x(2:obj.ParticleArray.Nt), iv_x(2:obj.ParticleArray.Nt)];
            jj_B = [ia_x(1:obj.ParticleArray.Nt-1), ia_x(1:obj.ParticleArray.Nt-1), ia_x(2:obj.ParticleArray.Nt)];
            vv_B = [0.5*obj.ParticleArray.dt*obj.ParticleArray.dt*unos(obj.ParticleArray.Nt-1), 0.5*obj.ParticleArray.dt*unos(obj.ParticleArray.Nt-1), 0.5*obj.ParticleArray.dt*unos(obj.ParticleArray.Nt-1)];

            numRows = 6*obj.ParticleArray.Nt;
            numCols = 3*obj.ParticleArray.Nt;
            B = sparse([ii_B, ii_B+obj.ParticleArray.Nt, ii_B+2*obj.ParticleArray.Nt],...
                [jj_B, jj_B+obj.ParticleArray.Nt, jj_B+2*obj.ParticleArray.Nt], ...
                [vv_B, vv_B, vv_B], numRows, numCols);
            
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
            
            B2 = sparse(idx1_vec, idx2_vec, val_vec, ii*numRows, ii*numCols);

            
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

            A0 = sparse([ii_A0, ii_A0+obj.ParticleArray.Nt, ii_A0+2*obj.ParticleArray.Nt], ...
                [jj_A0, jj_A0+1, jj_A0+2], ...
                [vv_A0, vv_A0, vv_A0], numRows, numCols);
            
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

            A02 = sparse(idx1_vec, idx2_vec, val_vec, ii*numRows, ii*numCols);


        end
        

        function [obj, S_p] = get_PrimalS(obj, Ix, Iy, Iz, systemMatrix)
%% Calculation of the Primal Matrix
% This function calculates the matrix of the Primal System (for a
% derivation with respect to the electric field).
                obj.D_I_1 = {};
                obj.D_I_2 = {};
                obj.D_I_3 = {};
                
                obj.Bp = {};
                
                for ii = 1:length(Ix)
                    obj.D_I_1{ii}          = [Ix{ii}*obj.Ex(:); Ix{ii}*obj.Er(:); Ix{ii}*obj.Ez(:)]; % [kg*m*s^-2*A^-1]
                    obj.D_I_2{ii}          = [Iy{ii}*obj.Ex(:); Iy{ii}*obj.Er(:); Iy{ii}*obj.Ez(:)];
                    obj.D_I_3{ii}          = [Iz{ii}*obj.Ex(:); Iz{ii}*obj.Er(:); Iz{ii}*obj.Ez(:)];
                    
                

                    grad_E_1 = zeros(1,3*obj.ParticleArray.Nt);
                    grad_E_2 = zeros(1,3*obj.ParticleArray.Nt);
                    grad_E_3 = zeros(1,3*obj.ParticleArray.Nt);

                    obj.Bp{ii} = B_matrix_primal(obj.ParticleArray.t_vec,grad_E_1, grad_E_2, grad_E_3, obj.D_I_1{ii}, obj.D_I_2{ii}, obj.D_I_3{ii});


                    S_p{ii} = systemMatrix + obj.Bp{ii};
                
                end
        end
        
        
        function [obj, S_p] = get_PrimalS2(obj, D_I_1, D_I_2, D_I_3, systemMatrix)
%% Calculation of the Primal Matrix
% This function calculates the matrix of the Primal System (for a
% derivation with respect to the electric field).

                %for ii = 1:length(Ix)
%                 Ex = repmat(obj.Ex(:),obj.ParticleArray.N_particles);
%                 Ey = repmat(obj.Ey(:),obj.ParticleArray.N_particles);
%                 Ez = repmat(obj.Ez(:),obj.ParticleArray.N_particles);
%                 
%                 D_I_1          = [Ix*Ex(:); Ix*Er(:); Ix*Ez(:)]; % [kg*m*s^-2*A^-1]
%                 D_I_2          = [Iy*Ex(:); Iy*Er(:); Iy*Ez(:)];
%                 D_I_3          = [Iz*Ex(:); Iz*Er(:); Iz*Ez(:)];
%                     
                
% 
%                grad_E_1 = zeros(1,3*obj.ParticleArray.Nt);
%                grad_E_2 = zeros(1,3*obj.ParticleArray.Nt);
%                grad_E_3 = zeros(1,3*obj.ParticleArray.Nt);

               obj.Bp2 = B_matrix_primal_Array(obj.ParticleArray.t_vec, D_I_1, D_I_2, D_I_3, obj.ParticleArray.N_particles);


               S_p = systemMatrix + obj.Bp2;
                
                %end
        end


        function [dGdEx, dGdEy, dGdEz, dGda] = getdGdE(obj,ii, accelMatrix, ...
            ix_x, ix_y, ix_z, accelInterpMatrix)
%% dGdEx, dGdEy, dGdEz Calculation
% This function calculates the derivatives of the objective function
% towards the components of the electric field at each grid point

            xv_dual = obj.ParticleArray.xv_dual(:,ii);
            dGda = xv_dual' * (-accelMatrix);

            dGdax = dGda(ix_x);
            dGdEx = reshape(dGdax * accelInterpMatrix{ii}, obj.resolution);

            dGday = dGda(ix_y);
            dGdEy = reshape(dGday * accelInterpMatrix{ii}, obj.resolution);

            dGdaz = dGda(ix_z);
            dGdEz = reshape(dGdaz * accelInterpMatrix{ii}, obj.resolution);
            
            %dGda = xv_dual2' * (-accelMatrix2)


        end
        
        function [obj] = getdGdE2(obj, accelMatrix, ...
            accelInterpMatrix, xv_dual)
%% dGdEx, dGdEy, dGdEz Calculation
% This function calculates the derivatives of the objective function
% towards the components of the electric field at each grid point

            %xv_dual = obj.ParticleArray.xv_dual(:,ii);
            obj.dGda = xv_dual' * (-accelMatrix);
            
            [ix_x, ix_y, ix_z, ~] = obj.ParticleArray.get_Index3D();
            
            mult = 1:obj.ParticleArray.N_particles;
            mult = mult * 3* obj.ParticleArray.Nt;
            mult = [0, mult(1:(end-1))];
            idx_x = [];
            idx_y = [];
            idx_z = [];
            
            for ii = 1:length(mult)
                
                idx_x = [idx_x, ix_x + mult(ii)];
                idx_y = [idx_y, ix_y + mult(ii)];
                idx_z = [idx_z, ix_z + mult(ii)];
                
            end
            
            obj.dGdEx2 = obj.dGda(idx_x) * accelInterpMatrix;
            obj.dGdEy2 = obj.dGda(idx_y) * accelInterpMatrix;
            obj.dGdEz2 = obj.dGda(idx_z) * accelInterpMatrix;
            
            obj.dGdEx = sum(reshape(obj.dGdEx2, prod(obj.resolution), []), 2);
            obj.dGdEy = sum(reshape(obj.dGdEy2, prod(obj.resolution), []), 2);
            obj.dGdEz = sum(reshape(obj.dGdEz2, prod(obj.resolution), []), 2);
            
            
            obj.dGdEx = reshape(obj.dGdEx, obj.resolution);
            obj.dGdEy = reshape(obj.dGdEy, obj.resolution);
            obj.dGdEz = reshape(obj.dGdEz, obj.resolution);

%             dGdax = dGda(ix_x);
%             dGdEx = reshape(dGdax * accelInterpMatrix, obj.resolution);
% 
%             dGday = dGda(ix_y);
%             dGdEy = reshape(dGday * accelInterpMatrix, obj.resolution);
% 
%             dGdaz = dGda(ix_z);
%             dGdEz = reshape(dGdaz * accelInterpMatrix, obj.resolution);
            
            %dGda = xv_dual2' * (-accelMatrix2)


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
