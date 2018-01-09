function [ energy, MMF_total, gap ] = energy_function( alpha_e, alpha_m, theta, gap, Coils, Winding, I, res )
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
    alpha_m = deg2rad(alpha_m);
    
    while alpha_m > 2*pi
        alpha_m = alpha_m - 2*pi;
    end
    
    while alpha_m < 0
        alpha_m = alpha_m + 2*pi;
    end
    
    [dif,shift] = min(abs(theta-alpha_m));
    gap = circshift(gap, shift);
    
    gap = gap.*10^-3; %meters
    
    gap_inv = gap.^-1;
    
    %% Handle electrical shift
    alpha_e = deg2rad(alpha_e);
    
    edges = 0:0.2618:6.2832;
    deg2coil = discretize(theta,edges);
    
    while alpha_e > 2*pi
        alpha_e = alpha_e - 2*pi;
    end

    while alpha_e < 0
        alpha_e = alpha_e + 2*pi;
    end

    [dif,index] = min(abs(theta-alpha_e));
    I_shift = deg2coil(index);
    I = circshift(I, I_shift);   
    
    MMF_total = MMF(theta, Coils, Winding, I);
    
    %% Calculate energy
    energy = ((Machine_height_m * mu_0) / 4) *...
        (trapz(((2*Stator_radius_m .*((MMF_total.^2).*gap_inv)) +...
        (MMF_total.^2)),2).*res);
    
end

