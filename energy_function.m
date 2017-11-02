function [ energy ] = energy_function( alpha1, alpha2, theta, gap, MMF_total, res )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   alpha1 -> electrical shift
%   alpha2 -> mechanical shift
    
    Machine_height = 2.5; %inches
    Machine_height_m = Machine_height*0.0254; %meters
    
    Stator_radius = 2.618; %inches
    Stator_radius_m = Stator_radius*0.0254; %meters
    
    mu_0 = 4*pi*10^-7; %henries/meter
    
    %% Handle mechanical shift
    alpha2 = deg2rad(alpha2);
    
    while alpha2 > 2*pi
        alpha2 = alpha2 - 2*pi;
    end
    
    while alpha2 < 0
        alpha2 = alpha2 + 2*pi;
    end
    
    [dif,shift] = min(abs(theta-alpha2));
    gap = circshift(gap, shift);
    
    gap = gap.*10^-3; %meters
    
    gap_inv = gap.^-1;
    
    %% Handle electrical shift
    alpha1 = deg2rad(alpha1);
    
    while alpha1 > 2*pi
        alpha1 = alpha1 - 2*pi;
    end

    while alpha1 < 0
        alpha1 = alpha1 + 2*pi;
    end

    [dif,shift] = min(abs(theta-alpha1));
    MMF_total = circshift(MMF_total, shift);   
    
    %% Calculate energy
    energy = ((Machine_height_m * mu_0) / 4) *...
        (trapz(((2*Stator_radius_m .*((MMF_total.^2).*gap_inv)) +...
        (MMF_total.^2)),2).*res);
end

