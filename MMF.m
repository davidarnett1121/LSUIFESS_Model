function [MMF_total] = MMF(theta, Coils, Winding, I)
%MMF Calculate MMF from each coil
%   Detailed explanation goes here

MMF = zeros(Coils,length(theta));

MMF = I.'.*Winding;

MMF_total = sum(MMF);
  
end

