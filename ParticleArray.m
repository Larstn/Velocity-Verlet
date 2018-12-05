classdef ParticleArray
    
    properties
        N_particles
        
        E_init
        v_init_abs
             
        nmass
        ncharges
        
        xyz_init
        v_init
        
        accelerationFunction
        accelerationFunction2
        
        Nt
        t_start
        t_end
        
        xv
        xv_dual
        
        accelerations

    end 
    
    properties (SetAccess=private, GetAccess=public, Hidden)
        t_vec 
        dt
    end
    
    properties (Dependent, Hidden)
                
        E_vec
        v_vec
        
        xx
        yy
        zz
        
        dxxdp
        dyydp
        dzzdp
        
        vx
        vy
        vz 
        
        dvxdp
        dvydp
        dvzdp
           
    end 
       
    properties (Constant)
       e = 1.60217662e-19;
       u = 1.6605e-27;
    end
    
    methods
        
        function obj = ParticleArray(Ev,Evstr,nm,nc,xyzi, tstart, tend, Nt, Nr)
            
            if nargin > 0
            %% Inputs: Energy/Velocity, String: 'Energy' or 'Velocity', n_charges, n_masses, initial positions, start time, end time, time steps
                assert(isa(nc,'double'),'Number of elementary charges must be a double')     
                assert(isa(nm,'double'),'Number of unit masses must be a double') 
                assert((isa(xyzi,'double') && isequal(size(xyzi),[3 , Nr])), 'Initial Position must be an row vector of doubles length 3')
                assert(strcmp(Evstr,'Energy') || strcmp(Evstr,'Velocity'), 'Allowed modes: Velocity or Energy')
                assert(isa(tstart,'double'),'Start Time must be a double') 
                assert(isa(tend,'double'),'End Time must be a double') 
                assert(isa(Nt,'double'),'Number of Timesteps must be a double') 
                
                obj.N_particles = Nr;

                if strcmp(Evstr,'Energy')
                     assert((isa(Ev,'double') && isequal(size(Ev),[3 ,Nr])), 'Energy Vector must contain two angles (size = 1,3)')

                     obj.E_init = Ev(1,:);
                     obj.v_init_abs = sqrt(2*Ev/(nm*u));

                     obj.v_init = zeros(Nr,3);
                     obj.v_init(1,:) = obj.v_init_abs*cos(Ev(2,:))*sin(Ev(3,:));
                     obj.v_init(2,:) = obj.v_init_abs*sin(Ev(2,:))*sin(Ev(3,:));
                     obj.v_init(3,:) = obj.v_init_abs*cos(Ev(3,:));

                elseif strcmp(Evstr,'Velocity')
                     assert((isa(Ev,'double') && isequal(size(Ev),[3 ,Nr])), 'Velocity Vector must contain 3 (carthesian) elements')

                     obj.v_init_abs = vecnorm(Ev,1);
                     obj.E_init = 0.5.*nm.*obj.u.*vecnorm(Ev,1).^2;

                     obj.v_init = zeros(3,Nr);
                     obj.v_init(1,:) = Ev(1,:);
                     obj.v_init(2,:) = Ev(2,:);
                     obj.v_init(3,:) = Ev(3,:);
                end

                obj.nmass = nm(1);
                obj.ncharges = nc(1);
                obj.xyz_init = xyzi;
                obj.t_start = tstart(1);
                obj.t_end = tend(1);
                obj.Nt = Nt(1);
            end
            
        end
        
        function [i_x, i_y, i_z, i_vx, i_vy, i_vz] = get_Index3D(obj)
            
            i_x = 1:obj.Nt(1);
            i_y = i_x+obj.Nt(1);
            i_z = i_y+obj.Nt(1);
            
            i_vx = i_z+obj.Nt(1);
            i_vy = i_vx+obj.Nt(1);
            i_vz = i_vy+obj.Nt(1);
            
        end
        
        function absv = get.v_vec(obj)
            
            [i_x, i_y, i_z, i_vx, i_vy, i_vz] = get_Index3D(obj);
            absv = sqrt(obj.xv(i_vx).^2 + obj.xv(i_vy).^2 + obj.xv(i_vz).^2);
            
        end
        
        function Evec = get.E_vec(obj)
            
            Evec = 0.5*obj.nmass*obj.v_vec.^2;
            
        end
        
        function xxx = get.xx(obj)
            
            [i_x, i_y, i_z, i_vx, i_vy, i_vz] = get_Index3D(obj);
            xxx = obj.xv(i_x,:);
            
        end
        
        function yyy = get.yy(obj)
            
            [i_x, i_y, i_z, i_vx, i_vy, i_vz] = get_Index3D(obj);
            yyy = obj.xv(i_y,:);
            
        end 
        
        function zzz = get.zz(obj)
            
            [i_x, i_y, i_z, i_vx, i_vy, i_vz] = get_Index3D(obj);
            zzz = obj.xv(i_z,:);
            
        end 
        
        function dxxxdp = get.dxxdp(obj)
            
            [i_x, i_y, i_z, i_vx, i_vy, i_vz] = get_Index3D(obj);
            dxxxdp = obj.xv_dual(i_x,:);
            
        end 
        
        function dyydp = get.dyydp(obj)
            
            [i_x, i_y, i_z, i_vx, i_vy, i_vz] = get_Index3D(obj);
            dyyydp = obj.xv_dual(i_y,:);
            
        end         
        
        function dzzzdp = get.dzzdp(obj)
            
            [i_x, i_y, i_z, i_vx, i_vy, i_vz] = get_Index3D(obj);
            dzzzdp = obj.xv_dual(i_z,:);
            
        end
        
        function vvx = get.vx(obj)
            
            [i_x, i_y, i_z, i_vx, i_vy, i_vz] = get_Index3D(obj);
            vvx = obj.xv(i_vx,:);
            
        end
        
        function vvy = get.vy(obj)
            
            [i_x, i_y, i_z, i_vx, i_vy, i_vz] = get_Index3D(obj);
            vvy = obj.xv(i_vy,:);
            
        end 
        
        function vvz = get.vz(obj)
            
            [i_x, i_y, i_z, i_vx, i_vy, i_vz] = get_Index3D(obj);
            vvz = obj.xv(i_vz,:);
            
        end 
        
        function dvvxdp = get.dvxdp(obj)
            
            [i_x, i_y, i_z, i_vx, i_vy, i_vz] = get_Index3D(obj);
            dvvxdp = obj.xv_dual(i_vx,:);
            
        end 
        
        function dvvydp = get.dvydp(obj)
            
            [i_x, i_y, i_z, i_vx, i_vy, i_vz] = get_Index3D(obj);
            dvvydp = obj.xv_dual(i_vy,:);
            
        end         
        
        function dvvzdp = get.dvzdp(obj)
            
            [i_x, i_y, i_z, i_vx, i_vy, i_vz] = get_Index3D(obj);
            dvvzdp = obj.xv_dual(i_vz,:);
            
        end  
        
        function tvec = get.t_vec(obj)
            
            %tvec = zeros(obj.Nt(1), obj.N_particles);
            %for ii = 1:obj.N_particles
                tvec = linspace(obj.t_start,obj.t_end,obj.Nt);
            %end
            
        end
        
        function dt = get.dt(obj)
            
            dt = (obj.t_end - obj.t_start)./(obj.Nt(1)-1);
            
        end  
        
        function obj = set.accelerationFunction(obj,fieldAccel)
            obj.accelerationFunction =  @(t, xyz) fieldAccel(t, xyz).*((obj.ncharges.*obj.e)/(obj.nmass.*obj.u));
        end
        
%         function obj = set.accelerationFunction2(obj,fieldAccel)
%             obj.accelerationFunction2 =  @(t, xyz) fieldAccel(t, xyz).*((obj.ncharges(1).*obj.e)/(obj.nmass(1).*obj.u));
%         end
        
        function obj = runTrajectory(obj)
            
            assert(~isempty(obj.accelerationFunction),...
                'Cannot run without an acceleration Function')
            xs = zeros(obj.Nt(1), 3, obj.N_particles);
            vs = zeros(obj.Nt(1), 3, obj.N_particles);
            xs(1,1:3,:) = obj.xyz_init; % m/s
            vs(1,1:3,:) = obj.v_init;   % m/s

            obj.accelerations = zeros(obj.Nt(1), 3, obj.N_particles);
            %currentAcceleration2 = zeros(150,3,5);
            currentAcceleration = obj.accelerationFunction(obj.t_vec(1,:),...
                permute(xs(1,1:3,:),[2,3,1]));
            
%             currentAcceleration2(1,:,1) = obj.accelerationFunction(obj.t_vec(1,:),...
%                 xs(1,1:3,1));
%             
%             currentAcceleration2(1,:,2) = obj.accelerationFunction(obj.t_vec(1,:),...
%                 xs(1,1:3,2));
%             
%             currentAcceleration2(1,:,3) = obj.accelerationFunction(obj.t_vec(1,:),...
%                 xs(1,1:3,3));
%             
%              currentAcceleration2(1,:,4) = obj.accelerationFunction(obj.t_vec(1,:),...
%                 xs(1,1:3,4));
%             
%             currentAcceleration2(1,:,5) = obj.accelerationFunction(obj.t_vec(1,:),...
%                 xs(1,1:3,5));
            
            
%             currentAcceleration2 = zeros(150,3,5);
%             
%             currentAcceleration2(1,:,1) = obj.accelerationFunction(1,...
%                 permute(xs(1,1:3,1),[2 1 3]));
%             
%             currentAcceleration2(1,:,2) = obj.accelerationFunction2(1,...
%                 permute(xs(1,1:3,2),[2 1 3]));
%             
%             currentAcceleration2(1,:,3) = obj.accelerationFunction2(1,...
%                 permute(xs(1,1:3,3),[2 1 3]));
%             
%             currentAcceleration2(1,:,4) = obj.accelerationFunction2(1,...
%                 permute(xs(1,1:3,4),[2 1 3]));
%             
%             currentAcceleration2(1,:,5) = obj.accelerationFunction2(1,...
%                 permute(xs(1,1:3,5),[2 1 3]));
           % assert(isequal(size(currentAcceleration), [3, 1]));

            obj.accelerations(1,1:3,:) = currentAcceleration;

            for nn = 1:obj.Nt(1)-1

                xs(nn+1,:,:) = xs(nn,:,:) + permute(permute(vs(nn,:,:),[2 3 1]).*obj.dt,[3 1 2]) + ...
                    0.5.*permute(currentAcceleration.*obj.dt.*obj.dt,[3 1 2]); % m
                nextAcceleration = obj.accelerationFunction(...
                    obj.t_vec(nn+1), permute(xs(nn+1,1:3,:), [2,3,1])); % m/s^2
                vs(nn+1,:,:) = vs(nn,:,:) + 0.5.*permute(obj.dt.*(currentAcceleration...
                    + nextAcceleration),[3 1 2]); % m/s
                obj.accelerations(nn+1,1:3,:) = nextAcceleration;

                currentAcceleration = nextAcceleration; %m?s^2
                
                
%                 currentAcceleration2(nn+1,:,1) = obj.accelerationFunction(obj.t_vec(1,:),...
%                 xs(nn+1,1:3,1));
%             
%             
%                 currentAcceleration2(nn+1,:,2) = obj.accelerationFunction(obj.t_vec(1,:),...
%                 xs(nn+1,1:3,2));
%             
%             
%                 currentAcceleration2(nn+1,:,3) = obj.accelerationFunction(obj.t_vec(1,:),...
%                 xs(nn+1,1:3,3));
%             
%             
%                 currentAcceleration2(nn+1,:,4) = obj.accelerationFunction(obj.t_vec(1,:),...
%                 xs(nn+1,1:3,4));
%             
%             
%                 currentAcceleration2(nn+1,:,5) = obj.accelerationFunction(obj.t_vec(1,:),...
%                 xs(nn+1,1:3,5));
            
            end

            obj.xv = permute([xs(:,1,:); xs(:,2,:); xs(:,3,:);
                vs(:,1,:); vs(:,2,:); vs(:,3,:)],[1 3 2]); %m 
            
        end 
         
        
        
    end
    
end 