function [ Sw ] = Sw( p,h ,coeff)

if ~exist('coeff','var')
    coeff = [1.00207; 4.42607; -5.47456; 5.02875; -1.24791; ...
    -2.26162; 4.38441; -1.79088; 3.69276; 5.17644;...
    7.30984; 1.29239; -1.00333; 3.9881; -9.09697; 1.29267; -6.28359; ...
    2.82282; -3.91952; 2.54342; -9.38879];
end
%rho_w(5), rho_s(5), h_w(7), h_s(4)
if exist('coeff','var')
    [ rhow_coeff,rhos_coeff,hw_coeff ,hs_coeff] = coeff_to_fun_coeff( coeff );
    if (h < hw_p(p,hw_coeff))
        Sw = 1;
    elseif (h > hs_p(p,hs_coeff))
        Sw = 0;
    else
        Sw = (rhos_p(p,rhos_coeff,hs_coeff)*(hs_p(p,hs_coeff)-h))/(h*(rhow_p(p,rhow_coeff,hw_coeff)-rhos_p(p,rhos_coeff,hs_coeff))...
            -(hw_p(p,hw_coeff)*rhow_p(p,rhow_coeff,hw_coeff)-hs_p(p,hs_coeff)*rhos_p(p,rhos_coeff,hs_coeff)));
    end
else
    if (h < hw_p(p))
        Sw = 1;
    elseif (h > hs_p(p))
        Sw = 0;
    else
        Sw = (rhos_p(p)*(hs_p(p)-h))/(h*(rhow_p(p)-rhos_p(p))-(hw_p(p)*rhow_p(p)-hs_p(p)*rhos_p(p)));
    end
end

end

