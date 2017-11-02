clear all
close all
clc

sample_f=1*10^6;
res = 1/sample_f;

mu_0 = 4*pi*10^-7; %henries/meter

theta = 0:res:2*pi; %radians


%% Define the turns function as a matrix for faster computation.

Coils = 24;
Coil_turns = 55;
Beta = pi/12;

turns = zeros(Coils,length(theta));

for x = 1:Coils
    if (x-1) <= 15
        for i = 1:length(theta)
            if (theta(i) >= (Beta/2 + 2*Beta + (x-1)*Beta)) &&...
                    (theta(i) < (Beta/2 + 8*Beta + (x-1)*Beta))
                turns(x,i) = Coil_turns;
            else
                turns(x,i) = 0;
            end
        end
    elseif (x-1) >= 22
        for i = 1:length(theta)
            if (theta(i) >= (Beta/2 + ((x-1)-22)*Beta)) &&...
                    (theta(i) < (Beta/2 + 6*Beta + ((x-1)-22)*Beta))
                turns(x,i) = Coil_turns;
            else
                turns(x,i) = 0;
            end
        end   
    else
        for i = 1:length(theta)
            if (theta(i) > (Beta/2 + ((x-1)-16)*Beta)) &&...
                    (theta(i) < (Beta/2 + 18*Beta + ((x-1)-16)*Beta))
                turns(x,i) = 0;
            else
                turns(x,i) = Coil_turns;
            end
        end     
    end
    
    if (x-1) == 6 || (x-1) == 7 || (x-1) == 8 || (x-1) == 9 ||...
            (x-1) == 10 || (x-1) == 11 || (x-1) == 18 || (x-1) == 19 ||...
            (x-1) == 20 || (x-1) == 21 || (x-1) == 22 || (x-1) == 23
        turns(x,:) = turns(x,:)*-1;
    end
end

% Create figure
figure('Name','Graphical Representation of the Turns Function')

    % Create plot
    plot(theta,turns(1,:),theta,turns(17,:));

    % Create xlabel
    xlabel('Electrical Angle (rad)');

    % Create ylabel
    ylabel('Turns');

    % Create x-limits of the axes
    xlim([0 6.2832]);

    
%% Calculate the average of the turns functions by each coil.

turns_avg = (1/(2*pi))*trapz(turns,2).*res; %


%% Define airgap 

gap_1 = gap_function(theta); %milimeters

% Solve for airgap translation
 
 x = 0;
 y = 0;
 
radial_shift = sqrt(x^2 + y^2);
shift_angle = atan2(y,x);
            
gap = gap_1 + radial_shift*sin(theta+shift_angle);

gap_avg = trapz(gap,2).*res*1E-3; %meter

% Create figure
figure('Name','Graphical Representation of the Air Gap')

    % Create plot
    plot(theta,gap);

    % Create xlabel
    xlabel('Electrical Angle (rad)');

    % Create ylabel
    ylabel('gap (mm)');

    % Create x-limits of the axes
    xlim([0 6.2832]);
    
gap_inv = gap.^-1;

gap_inv_avg = (1/(2*pi))*trapz(gap_inv,2); 


%% Define the modified winding fuction as a matrix.

 M_Winding = zeros(Coils,length(theta));
 
 for x = 1:Coils
     M_Winding(x,:) = turns(x,:)-(1/(2*pi*gap_inv_avg))*(trapz(gap_inv.*turns(x,:),2));
 end


%% Coil Order
% 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 
% Q Q D D D D Q Q D D  D  D  Q  Q  D  D  D  D  Q  Q  D  D  D  D

I_d = 0.42; %Amps
I_q = 0.952*I_d; %Amps

I = [I_q I_q I_d I_d I_d I_d I_q I_q I_d I_d I_d I_d I_q I_q I_d I_d...
                                          I_d I_d I_q I_q I_d I_d I_d I_d];

%% Calculate MMF from each coil.

MMF = zeros(Coils,length(theta));

MMF = I.'.*M_Winding;

MMF_total = sum(MMF);

% Create figure
figure('Name','Magneto-Motive Force')

    % Create plot
    plot(theta,MMF_total);

    % Create xlabel
    xlabel('Electrical Angle (theta)');

    % Create ylabel
    ylabel('Modified MMF (A)');

    % Create limits of the axes
    xlim([0 6.2832]);
    ylim([-1000 1000]);

  
%% Solve for energy in the airgap

energy = zeros(1,361);

for alpha = 0:360
    energy(alpha+1) = energy_function( alpha, 0, theta, gap, MMF_total, res );
end

alpha = 0:360;

figure('Name','Air Gap Energy')

    % Create plot
    plot(alpha,energy);

    % Create xlabel
    xlabel('Mechanical Rotor Angle (deg)');

    % Create ylabel
    ylabel('Energy (J)');

    % Create x-limits of the axes
    xlim([0 360]);
    
    
%%Solve for torque

torque = diff(energy)/deg2rad(1);

torque = [0 torque];

figure('Name','Torque')

    % Create plot
    plot(alpha,torque);

    % Create xlabel
    xlabel('Mechanical Rotor Angle (deg)');

    % Create ylabel
    ylabel('Torque (N*m)');

    % Create x-limits of the axes
    xlim([0 360]);
    

%% Solve for translational force

 x = 0;
 y = 0;
 
    radial_shift = sqrt(x^2 + y^2);
    shift_angle = atan2(y,x);
    
    alpha = deg2rad(0);
    
    while alpha > 2*pi
        alpha = alpha - 2*pi;
    end
    
    while alpha < 0
        alpha = alpha + 2*pi;
    end
    
    [dif,shift] = min(abs(theta-alpha));
    gap = circshift(gap_1, shift);

    gap = gap.*1E-3 + radial_shift*sin(theta+shift_angle).*1E-3;

    % Create figure
    figure('Name','Graphical Representation of the Air Gap')

        % Create plot
        plot(theta,gap);

        % Create xlabel
        xlabel('Electrical Angle (rad)');

        % Create ylabel
        ylabel('gap (m)');

        % Create x-limits of the axes
        xlim([0 6.2832]);

%%
    
    gap_avg = trapz(gap,2).*res; %meter

    gap_inv = gap.^-1;

    gap_inv_avg = (1/(2*pi))*trapz(gap_inv,2); 
    
     M_Winding = zeros(Coils,length(theta));
 
     for x = 1:Coils
         M_Winding(x,:) = turns(x,:)-(1/(2*pi*gap_inv_avg))*(trapz(gap_inv.*turns(x,:),2));
     end
 
    MMF = I.'.*M_Winding;

    MMF_total = sum(MMF);
    
    % Solve for energy in the airgap
    
    H = MMF_total.*gap_inv;
    
    B = mu_0.*H;
    
%     F_x = zeros(1,361);
%     F_y = zeros(1,361);
%     
%     for alpha = 0:360
%         [dif,loc] = min(abs(theta-deg2rad(alpha)));
%         F_x(alpha+1) = (((B(loc)^2)*66.4718*49.784)/(2*mu_0)).*cos(theta(loc)+deg2rad(alpha));
%         F_y(alpha+1) = (((B(loc)^2)*66.4718*49.784)/(2*mu_0)).*sin(theta(loc)+deg2rad(alpha));
%     end
% 
     alpha = 0:360;
    
%     th = deg2rad(45);
%     
     F_x = (((B.^2).*0.034804553091815.*0.0635)./(2*mu_0)).*cos(theta + deg2rad(45));

     Fx_total = trapz(gap,F_x,2);

     F_y = (((B.^2).*0.034804553091815.*0.0635)./(2*mu_0)).*sin(theta + deg2rad(45));

     Fy_total = trapz(gap,F_y,2);
    
 figure('Name','Force')
 
     % Create plot
     plot(theta, F_x, theta, F_y,...
          theta, gap.*10^4.26,...
          theta, MMF_total.*10^-0.28,...
          theta, B.*10^2.62);
      
     legend('F_x', 'F_y', 'gap', 'MMF', 'B');

     % Create xlabel
     xlabel('Mechanical Rotor Angle (rad)');
 
     % Create ylabel
     ylabel('Force (N)');
 
     % Create x-limits of the axes
     xlim([0 2*pi]);
     
    