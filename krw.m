function [ krw ] = krw( p,h ,coeff)


% Slstar = (Sw(p,h)-0.3)/(1-0.3-0.05);
% krw = Slstar^4;
krw = Sw(p,h,coeff);
end

