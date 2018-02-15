function [ rhow ] = rhow_p( p,rho_coeff,h_coeff )

rhoconv = 1E3;
hconv = 1E4;
pconv = 10;
if exist('rho_coeff','var')
    if exist('h_coeff','var')
        h = hw_p(p,h_coeff);
    else
        h = hw_p(p);
    end
    rhow = rho_coeff(1) + rho_coeff(2)*(p.*pconv) ...
                + rho_coeff(3)*(h.*hconv) ...
                + rho_coeff(4)*(p.*pconv.*h.*hconv) ...
                + rho_coeff(5)*(h.*hconv).^2;
else
    h = hw_p(p);
    rhow = 1.00207 + 4.42607e-11*(p.*pconv) ...
                - 5.47456e-12*(h.*hconv) ...
                + 5.02875e-21*(p.*pconv.*h.*hconv) ...
                -1.24791e-21*(h.*hconv).^2;
end

rhow = rhow.*rhoconv;
end

