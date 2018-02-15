function [ rhow_coeff,rhos_coeff,hw_coeff ,hs_coeff] = coeff_to_fun_coeff( coeff_in )

coeffPowers = [1; 1e-11; 1e-12; 1e-21; 1e-21; ...
    1e-5; 1e-9; 1e-19; 1e-36; 1e-41;...
    1e9; 1e2; 1e-6; 1e-15; 1e15; 1e22; 1e27; ...
    1e10; 1e15; 1e21; 1e-8];
coeff = coeff_in.*coeffPowers;



rhow_coeff = coeff(1:5);
rhos_coeff = coeff(6:10);
hw_coeff = coeff(11:17);
hs_coeff = coeff(18:21);

end

