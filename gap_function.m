function [ gap ] = gap_function( theta )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    d_mid = 155.03; %milimeters
    r_mid = d_mid/2; %milimeters
    r_inner = 134.92/2; %milimeters

    gap_interpole = r_mid - r_inner + 1; %milimeters
    gap_pole = 1; %milimeter
    
    rad_30 = deg2rad(30);
    rad_60 = deg2rad(60);
    rad_120 = deg2rad(120);
    rad_150 = deg2rad(150);
    rad_210 = deg2rad(210);
    rad_240 = deg2rad(240);
    rad_300 = deg2rad(300);
    rad_330 = deg2rad(330);

    gap = zeros(1,length(theta));

    for i = 1:length(theta)
        if (((theta(i) >= rad_30) && (theta(i) < rad_60)) ||...
            ((theta(i) >= rad_120) && (theta(i) < rad_150)) ||...
            ((theta(i) >= rad_210) && (theta(i) < rad_240)) ||...
            ((theta(i) >= rad_300) && (theta(i) < rad_330)))
        
            gap(i) = gap_pole;
        
        else
        
            gap(i) = gap_interpole;
        
        end
    end

end

