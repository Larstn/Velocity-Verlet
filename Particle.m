classdef Particle
    
    properties 
        E_init
        v_init_abs
             
        nmass
        ncharges
        
        xyz_init
        v_init
        
        accelerationFunction
        
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
        function obj = Particle(Ev,Evstr,nm,nc,xyzi, tstart, tend, Nt)
            
            if nargin > 0
            %% Inputs: Energy/Velocity, String: 'Energy' or 'Velocity', n_charges, n_masses, initial positions, start time, end time, time steps
                assert(isa(nc,'double'),'Number of elementary charges must be a double')     
                assert(isa(nm,'double'),'Number of unit masses must be a double') 
                assert((isa(xyzi,'double') && isequal(size(xyzi),[3 ,1])), 'Initial Position must be an row vector of doubles length 3')
                assert(strcmp(Evstr,'Energy') || strcmp(Evstr,'Velocity'), 'Allowed modes: Velocity or Energy')
                assert(isa(tstart,'double'),'Start Time must be a double') 
                assert(isa(tend,'double'),'End Time must be a double') 
                assert(isa(Nt,'double'),'Number of Timesteps must be a double') 


                if strcmp(Evstr,'Energy')
                     assert((isa(Ev,'double') && isequal(size(Ev),[3 ,1])), 'Energy Vector must contain two angles (size = 1,3)')

                     obj.E_init = Ev(1);
                     obj.v_init_abs = sqrt(2*Ev/(nm*u));

                     obj.v_init = zeros(1,3);
                     obj.v_init(1) = obj.v_init_abs*cos(Ev(2))*sin(Ev(3));
                     obj.v_init(2) = obj.v_init_abs*sin(Ev(2))*sin(Ev(3));
                     obj.v_init(3) = obj.v_init_abs*cos(Ev(3));

                elseif strcmp(Evstr,'Velocity')
                     assert((isa(Ev,'double') && isequal(size(Ev),[3 ,1])), 'Velocity Vector must contain 3 (carthesian) elements (size = 1,3)')

                     obj.v_init_abs = norm(Ev);
                     obj.E_init = 0.5*nm*obj.u*norm(Ev)^2;

                     obj.v_init = zeros(1,3);
                     obj.v_init(1) = Ev(1);
                     obj.v_init(2) = Ev(2);
                     obj.v_init(3) = Ev(3);
                end

                obj.nmass = nm;
                obj.ncharges = nc;
                obj.xyz_init = xyzi;
                obj.t_start = tstart;
                obj.t_end = tend;
                obj.Nt = Nt;
            end
            
        end
        
        function [i_x, i_y, i_z, i_vx, i_vy, i_vz] = get_Index3D(obj)
            
            i_x = 1:obj.Nt;
            i_y = i_x+obj.Nt;
            i_z = i_y+obj.Nt;
            
            i_vx = i_z+obj.Nt;
            i_vy = i_vx+obj.Nt;
            i_vz = i_vy+obj.Nt;
            
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
            xxx = obj.xv(i_x);
            
        end
        
        function yyy = get.yy(obj)
            
            [i_x, i_y, i_z, i_vx, i_vy, i_vz] = get_Index3D(obj);
            yyy = obj.xv(i_y);
            
        end 
        
        function zzz = get.zz(obj)
            
            [i_x, i_y, i_z, i_vx, i_vy, i_vz] = get_Index3D(obj);
            zzz = obj.xv(i_z);
            
        end 
        
        function dxxxdp = get.dxxdp(obj)
            
            [i_x, i_y, i_z, i_vx, i_vy, i_vz] = get_Index3D(obj);
            dxxxdp = obj.xv_dual(i_x);
            
        end 
        
        function dyydp = get.dyydp(obj)
            
            [i_x, i_y, i_z, i_vx, i_vy, i_vz] = get_Index3D(obj);
            dyyydp = obj.xv_dual(i_y);
            
        end         
        
        function dzzzdp = get.dzzdp(obj)
            
            [i_x, i_y, i_z, i_vx, i_vy, i_vz] = get_Index3D(obj);
            dzzzdp = obj.xv_dual(i_z);
            
        end
        
        function vvx = get.vx(obj)
            
            [i_x, i_y, i_z, i_vx, i_vy, i_vz] = get_Index3D(obj);
            vvx = obj.xv(i_vx);
            
        end
        
        function vvy = get.vy(obj)
            
            [i_x, i_y, i_z, i_vx, i_vy, i_vz] = get_Index3D(obj);
            vvy = obj.xv(i_vy);
            
        end 
        
        function vvz = get.vz(obj)
            
            [i_x, i_y, i_z, i_vx, i_vy, i_vz] = get_Index3D(obj);
            vvz = obj.xv(i_vz);
            
        end 
        
        function dvvxdp = get.dvxdp(obj)
            
            [i_x, i_y, i_z, i_vx, i_vy, i_vz] = get_Index3D(obj);
            dvvxdp = obj.xv_dual(i_vx);
            
        end 
        
        function dvvydp = get.dvydp(obj)
            
            [i_x, i_y, i_z, i_vx, i_vy, i_vz] = get_Index3D(obj);
            dvvydp = obj.xv_dual(i_vy);
            
        end         
        
        function dvvzdp = get.dvzdp(obj)
            
            [i_x, i_y, i_z, i_vx, i_vy, i_vz] = get_Index3D(obj);
            dvvzdp = obj.xv_dual(i_vz);
            
        end  
        
        function tvec = get.t_vec(obj)
            
            tvec = linspace(obj.t_start,obj.t_end,obj.Nt);
            
        end
        
        function dt = get.dt(obj)
            
            dt = (obj.t_end - obj.t_start)./(obj.Nt-1);
            
        end  
        
        function obj = set.accelerationFunction(obj,fieldAccel)
            obj.accelerationFunction =  @(t, xyz) fieldAccel(t, xyz)*((obj.ncharges*obj.e)/(obj.nmass*obj.u));
        end
        
        function obj = runTrajectory(obj)
            
            assert(~isempty(obj.accelerationFunction),...
                'Cannot run without an acceleration Function')
            xs = zeros(obj.Nt, 3);
            vs = zeros(obj.Nt, 3);
            xs(1,1:3) = obj.xyz_init; % m/s
            vs(1,1:3) = obj.v_init;   % m/s

            obj.accelerations = zeros(obj.Nt, 3);

            currentAcceleration = obj.accelerationFunction(obj.t_vec(1),...
                xs(1,1:3));
            assert(isequal(size(currentAcceleration), [3, 1]));

            obj.accelerations(1,1:3) = currentAcceleration;

            for nn = 1:obj.Nt-1

                xs(nn+1,:) = xs(nn,:) + vs(nn,:)*obj.dt + ...
                    0.5*obj.dt*obj.dt*currentAcceleration'; % m
                nextAcceleration = obj.accelerationFunction(...
                    obj.t_vec(nn+1), xs(nn+1,:)); % m/s^2
                vs(nn+1,:) = vs(nn,:) + 0.5*obj.dt*(currentAcceleration...
                    + nextAcceleration)'; % m/s
                obj.accelerations(nn+1,1:3) = nextAcceleration;

                currentAcceleration = nextAcceleration; %m?s^2
            end

            obj.xv = [xs(:,1); xs(:,2); xs(:,3);
                vs(:,1); vs(:,2); vs(:,3)]; %m 
            
         end 
    end
    
end 