classdef gapClass
    properties
        pole_gap
        interpole_gap
        theta
        array
        avg
        array_inv
        inv_avg
    end
    
    methods
        
        function obj = gapClass( a, b, res )
            
                d_mid = a; %milimeters
                r_mid = d_mid/2; %milimeters
                r_inner = b/2; %milimeters
                
                p1_s = deg2rad(30); %Starting angle of pole 1
                p1_e = deg2rad(60); %Ending angle of pole 1
                p2_s = deg2rad(120); %Starting angle of pole 2
                p2_e = deg2rad(150); %Ending angle of pole 2
                p3_s = deg2rad(210); %Starting angle of pole 3
                p3_e = deg2rad(240); %Ending angle of pole 3
                p4_s = deg2rad(330); %Starting angle of pole 4
                p4_e = deg2rad(360); %Ending angle of pole 4
                
                theta = 0:res:2*pi; %radians

                obj.interpole_gap = r_mid - r_inner + 1; %milimeters
                obj.pole_gap = 1; %milimeter

                obj.array = zeros(1,length(theta));

                for i = 1:length(theta)
                    if (((theta(i) >= p1_s) && (theta(i) < p1_e)) ||...
                        ((theta(i) >= p2_s) && (theta(i) < p2_e)) ||...
                        ((theta(i) >= p3_s) && (theta(i) < p3_e)) ||...
                        ((theta(i) >= p4_s) && (theta(i) < p4_e)))

                        obj.array(i) = obj.pole_gap;

                    else

                        obj.array(i) = obj.interpole_gap;

                    end
                    
                end
                
                obj.avg = (1/(2*pi))*trapz(obj.array,2)*res; 
                
                obj.array_inv = obj.array.^-1;
                
                obj.inv_avg = (1/(2*pi))*trapz(obj.array_inv,2)*res; 
                
        end
        
        function out = rotate_gap( obj, alpha )
            
                if alpha >= (2*pi)
                    alpha = alpha - (2*pi);
                end
                
                if alpha == 0
                    out = obj.array;
                end
                
                [dif,shift] = min(abs(obj.theta-alpha));
                out = circshift(obj.array, shift);
                
        end
        
        function out = translate_gap( obj, x, y )
            
            radial_shift = sqrt(x^2 + y^2);
            shift_angle = atan(y/x);
            
            gap = obj.array + radial_shift*sin(obj.theta+shift_angle);
            
            [gap_min, min_angle] = min(gap);
            
            if gap_min <= 0
                error('Collision Detected at %d degrees', rad2deg(min_angle))
            else
                out = gap;
            end        
                      
        end
        
        function out = mechanical_distortion( obj, omega )
              
            error('Mechanical deformation of the gap is not defined')
            
        end
        
        function out = saturation_effects( obj, i )
        
            error('Saturation effect on the gap is not defined')
            
        end
        
        function plot(obj)
            
            % Create figure
            figure('Name','Graphical Representation of the Air Gap')

                % Create plot
                plot(obj.theta,obj.array);

                % Create xlabel
                xlabel('Electrical Angle (rad)');

                % Create ylabel
                ylabel('gap (mm)');

                % Create x-limits of the axes
                xlim([0 (2*pi)]);
                
        end
        
    end
    
%   methods (Static)
%      function obj = createObj
%         prompt = {'Enter the Radius'};
%         dlgTitle = 'Radius';
%         rad = inputdlg(prompt,dlgTitle);
%         r = str2double(rad{:});
%         obj = CircleArea(r);
%      end
   end
    
    
        