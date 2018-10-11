function varargout = interplars(varargin)

    %% interplars

    % input: 1, 2 or 3 arrays, ascending order, equally spaced, value matrix
    % (1D,2D,3D) and then 1, 2 or 3 query arrays
    
    
    
    
    if nargin < 7
        if nargin == 3
            
            x_p = varargin{1};
            x_q = varargin{3};
            val = varargin{2};
            
            dx = x_p(2) - x_p(1);
            
            i_x = floor( (x_q - x_p(1)) ./ dx) + 1;
            
            %m = (val(2:end) - val(1:(end-1)))/dx;
            
            %b = val(1:(end-1)) - m.*x_p(1:(end-1));
            
            %m = m_vec;
            %b = b_vec;
            
            try
                varargout{1} = (val(i_x).*(x_p(i_x+1) - x_q) + val(i_x+1).*(x_q - x_p(i_x))) ./dx;
            catch 
                warning("Particles left the field")
                l = length(x_p);
                i_x( (i_x < 1) | (i_x > l) ) = l + 1;
                x_p = [x_p 0 0];
                val = [val 0 0];
                varargout{1} = (val(i_x).*(x_p(i_x+1) - x_q) + val(i_x+1).*(x_q - x_p(i_x))) ./dx;
            end
        else
            x_p = varargin{1};
            x_q = varargin{4};
            y_p = varargin{2};
            y_q = varargin{5};
            val = varargin{3};
            
            dx = x_p(2) - x_p(1);
            i_x = floor( (x_q - x_p(1)) ./ dx) + 1;
            
            dy = y_p(2) - y_p(1);
            i_y = floor( (y_q - y_p(1)) ./ dy) + 1;
            

            try
                fxy1 = ((x_p(i_x+1) - x_q).*val(i_x,i_y) + (x_q - x_p(i_x)).*val(i_x+1,i_y)) ./ dx;
                fxy2 = ((x_p(i_x+1) - x_q).*val(i_x,i_y+1) + (x_q - x_p(i_x)).*val(i_x+1,i_y+1)) ./dx;
                varargout{1} = ((y_p(i_y+1) - y_q).*fxy1 + (y_q - y_p(i_y)).*fxy2) /dy;
            catch e
                fprintf('%s \n',e.identifier)
                warning('2D out of box')
                x_p = [x_p 0 0];
                y_p = [y_p 0 0];
                val = padarray(val, [2 2], 0, 'post');
                
                l1 = length(x_p) - 2;
                l2 = length(y_p) - 2;
                
                i_oob = (i_x < 1) | (i_x > l1) | (i_y < 1) | (i_y > l2);
                i_x(i_oob) = l1 +1;
                i_y(i_oob) = l2 +1;
                
                
                fxy1 = (x_p(i_x+1) - x_q).*val(i_x,i_y) + (x_q - x_p(i_x)).*val(i_x+1,i_y);
                fxy2 = (x_p(i_x+1) - x_q).*val(i_x,i_y+1) + (x_q - x_p(i_x)).*val(i_x+1,i_y+1);
                varargout{1} = (y_p(i_y+1) - y_q).*fxy1 + (y_q - y_p(i_y)).*fxy2;
           % varargout{}
            
            end
        end
        else
            x_p = varargin{1};
            y_p = varargin{2};
            z_p = varargin{3};
            val = varargin{4};
            x_q = varargin{5};
            y_q = varargin{6};
            z_q = varargin{7};
            
            dx = x_p(2) - x_p(1);
            dy = y_p(2) - y_p(1);
            dz = z_p(2) - z_p(1);
            
            i_x = floor( (x_q - x_p(1)) ./ dx) + 1;
            i_y = floor( (y_q - y_p(1)) ./ dy) + 1;
            i_z = floor( (z_q - z_p(1)) ./ dz) + 1;
            
            try
                xd = (x_q - x_p(i_x)) ./dx;
                yd = (y_q - y_p(i_y)) ./dy;
                zd = (z_q - z_p(i_z)) ./dz;
                
                c00 = val(i_x,i_y,i_z)*(1-xd) + val(i_x+1,i_y,i_z)*xd;
                c01 = val(i_x,i_y,i_z+1)*(1-xd) + val(i_x+1,i_y,i_z+1)*xd;
                c10 = val(i_x,i_y+1,i_z)*(1-xd) + val(i_x+1,i_y+1,i_z)*xd;
                c11 = val(i_x,i_y+1,i_z+1)*(1-xd) + val(i_x+1,i_y+1,i_z+1)*xd;
                
                c0 = c00*(1-yd) + c10*yd;
                c1 = c01*(1-yd) + c11*yd;
                
                varargout{1} = c0*(1-zd) + c1*zd;

                
            catch
                warning('left field')
                x_p = [x_p 0 0];
                y_p = [y_p 0 0];
                z_p = [z_p 0 0];
                val = padarray(val, [2 2 2], 0, 'post');
                
                l1 = length(x_p) - 2;
                l2 = length(y_p) - 2;
                l3 = length(z_p) - 2;
                
                i_oob = (i_x < 1) | (i_x > l1) | (i_y < 1) | (i_y > l2) | (i_z < 1) | (i_z > 1);
                i_x(i_oob) = l1 +1;
                i_y(i_oob) = l2 +1;
                i_z(i_oob) = l3 +1;
                
                xd = (x_q - x_p(i_x)) ./dx;
                yd = (y_q - y_p(i_y)) ./dy;
                zd = (z_q - z_p(i_z)) ./dz;
                
                c00 = val(i_x,i_y,i_z)*(1-xd) + val(i_x+1,i_y,i_z)*xd;
                c01 = val(i_x,i_y,i_z+1)*(1-xd) + val(i_x+1,i_y,i_z+1)*xd;
                c10 = val(i_x,i_y+1,i_z)*(1-xd) + val(i_x+1,i_y+1,i_z)*xd;
                c11 = val(i_x,i_y+1,i_z+1)*(1-xd) + val(i_x+1,i_y+1,i_z+1)*xd;
                
                c0 = c00*(1-yd) + c10*yd;
                c1 = c01*(1-yd) + c11*yd;
                
                varargout{1} = c0*(1-zd) + c1*zd;
            end
            
            
    end
    % 1. Step: which queries are out of the boundaries? all those points, zero
    % 2. Step: if not all, figure out the linear function between all of the
    % array points by substracting and dividing 
    % 3. Step: decide in which box each query point belongs.. division instead
    % of sorting? Modulo 
    % 4. Step: Assign each query the actual value
    % 5. Step: Done

    % replace test inputs with real inputs 

end