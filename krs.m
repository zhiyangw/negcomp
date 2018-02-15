function [ krs ] = krs( p,h ,coeff)

% Slstar = (Sw(p,h)-0.3)/(1-0.3-0.05);
% krs = (1-Slstar)^2*(1-Slstar^2);
krs = 1-Sw(p,h,coeff);
end

