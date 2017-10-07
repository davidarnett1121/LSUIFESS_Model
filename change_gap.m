function [ gap ] = change_gap( gap, theta, alpha, x, y )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    [dif,shift] = min(abs(theta-alpha));
    gap = circshift(gap, shift);
    
    r = sqrt(x^2+y^2);
    alpha = atan(y/x);
    
    [dif,shift] = min(abs(theta-alpha));

end

