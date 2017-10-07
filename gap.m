function [ gap ] = gap( theta, alpha )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    d_mid = 155.03; %milimeters
    r_mid = d_mid/2; %milimeters
    r_inner = 134.92/2; %milimeters

    gap_interpole = r_mid - r_inner + 1; %milimeters
    gap_pole = 1; %milimeter

    gap = zeros(1,length(theta));

    for i = 1:length(theta)
        if (((theta(i) >= (deg2rad(30) + alpha)) && (theta(i) < (deg2rad(60) + alpha))) ||...
            ((theta(i) >= (deg2rad(120) + alpha)) && (theta(i) < (deg2rad(150) + alpha))) ||...
            ((theta(i) >= (deg2rad(210) + alpha)) && (theta(i) < (deg2rad(240) + alpha))) ||...
            ((theta(i) >= (deg2rad(330) + alpha)) && (theta(i) < (deg2rad(360) + alpha))))
        
            gap(i) = gap_pole;
        
        else
        
            gap(i) = gap_interpole;
        
        end
    end

end

