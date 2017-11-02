function [ gap ] = gap_function( theta )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

        d_mid = 155.03; %milimeters
        r_mid = d_mid/2; %milimeters
        r_inner = 134.92/2; %milimeters

        pole_angle = deg2rad(30);
        pole_offset = deg2rad(30);
        interpole_angle = deg2rad(60);

        p1_s = pole_offset; %Starting angle of pole 1
        p1_e = p1_s + pole_angle; %Ending angle of pole 1
        p2_s = p1_e + interpole_angle; %Starting angle of pole 2
        p2_e = p2_s + pole_angle; %Ending angle of pole 2
        p3_s = p2_e + interpole_angle; %Starting angle of pole 3
        p3_e = p3_s + pole_angle; %Ending angle of pole 3
        p4_s = p3_e + interpole_angle; %Starting angle of pole 4
        p4_e = p4_s + pole_angle; %Ending angle of pole 4

        interpole_gap = r_mid - r_inner + 1; %milimeters
        pole_gap = 1; %milimeter

        gap = zeros(1,length(theta));

        for i = 1:length(theta)
            if (((theta(i) >= p1_s) && (theta(i) < p1_e)) ||...
                ((theta(i) >= p2_s) && (theta(i) < p2_e)) ||...
                ((theta(i) >= p3_s) && (theta(i) < p3_e)) ||...
                ((theta(i) >= p4_s) && (theta(i) < p4_e)))

                gap(i) = pole_gap;

            else

                gap(i) = interpole_gap;

            end

        end

end

