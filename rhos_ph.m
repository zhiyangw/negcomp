function [rhos] = rhos_ph(p,h,rho_coeff,h_coeff ,usesat)

rhoconv = 1E3;
hconv = 1E4;
pconv = 10;
if ~exist('usesat','var')
    sat = Sw(p,h);
    if sat >= 1
        rhos = 0;
        return;
    elseif (sat>0 && sat < 1)
        if exist('h_coeff','var')
            h = hs_p(p,h_coeff);
        else
            h = hs_p(p);
        end
    end
end

if exist('rho_coeff','var')
    rhos = rho_coeff(1) + rho_coeff(2)*p.*pconv ...
        + rho_coeff(3).*p.*pconv.*h.*hconv ...
        + rho_coeff(4).*(p.*pconv).^4 ...
        + rho_coeff(5).*p.*pconv.*(h.*hconv).^3;
else
    rhos = -2.26162e-5 +4.38441e-9*p.*pconv ...
        - 1.79088e-19.*p.*pconv.*h.*hconv ...
        + 3.69276e-36.*(p.*pconv).^4 ...
        + 5.17644e-41.*p.*pconv.*(h.*hconv).^3;
end
rhos = rhos*rhoconv;


end
