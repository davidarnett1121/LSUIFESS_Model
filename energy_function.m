function [ energy ] = energy_function( alpha, theta, gap, MMF_total, res )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
    Machine_height = 2.5; %inches
    Machine_height_m = Machine_height*0.0254; %meters
    
    Stator_radius = 2.618; %inches
    Stator_radius_m = Stator_radius*0.0254; %meters
    
    mu_0 = 4*pi*10^-7; %henries/meter
        
    alpha = deg2rad(alpha);
    
    while alpha > 2*pi
        alpha = alpha - 2*pi;
    end
    
    while alpha < 0
        alpha = alpha + 2*pi;
    end
    
    [dif,shift] = min(abs(theta-alpha));
    gap = circshift(gap, shift);
    
    gap = gap.*10^-3; %meters
    
    gap_inv = gap.^-1;
    
    energy = ((Machine_height_m * mu_0) / 4) *...
        (trapz(((2*Stator_radius_m .*((MMF_total.^2).*gap_inv)) +...
        (MMF_total.^2)),2).*res);
end

