function [ viscs ] = viscs( p,h ,coeff)
if ~exist('coeff','var')
    
coeff = [1.00207; 4.42607; -5.47456; 5.02875; -1.24791; ...
    -2.26162; 4.38441; -1.79088; 3.69276; 5.17644;...
    7.30984; 1.29239; -1.00333; 3.9881; -9.09697; 1.29267; -6.28359; ...
    2.82282; -3.91952; 2.54342; -9.38879];
end
viscs = 0.1*1e-6*(0.407*(T_ph(p,h,coeff)-273.15)+80.4);

end

