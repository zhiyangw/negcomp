function [ pvars ] = fluidProperties(p,x2,prop )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if strcmpi(prop.primaryVariable,'pnv')
    error('wrong fluid properties')
end

 [PRES,ENTH,TEMP,SAT_W,SAT_S,RHO_W,RHO_S,RHO,H_W,H_S,VISC_W,VISC_S,KR_W,KR_S,NVAR] = offset();

pvars = zeros(NVAR,1);

prop.rhodata.coeff =[1.00207; 4.42607; -5.47456; 5.02875; -1.24791; ...
    -2.26162; 4.38441; -1.79088; 3.69276; 5.17644;...
    7.30984; 1.29239; -1.00333; 3.9881; -9.09697; 1.29267; -6.28359; ...
    2.82282; -3.91952; 2.54342; -9.38879];

[ rhow_coeff,rhos_coeff,hw_coeff ,hs_coeff] = coeff_to_fun_coeff( prop.rhodata.coeff );
coeff = prop.rhodata.coeff;

g = 9.81;


pvars(PRES) = pconv(p,prop);



if strcmpi(prop.primaryVariable, 'S');
    pvars(ENTH) = h_pS(pvars(PRES),x2,coeff);
else
    pvars(ENTH) = hconv(x2,prop);
end

pvars(SAT_W) = Sw(pvars(PRES),pvars(ENTH),coeff);
pvars(TEMP) = T_ph(pvars(PRES),pvars(ENTH),coeff);
pvars(VISC_W) = viscw(pvars(PRES),pvars(ENTH),coeff);
pvars(VISC_S) = viscs(pvars(PRES),pvars(ENTH),coeff);
pvars(KR_W) = krw(pvars(PRES),pvars(ENTH),coeff);
pvars(KR_S) = krs(pvars(PRES),pvars(ENTH),coeff);

pvars(RHO_W) = rhow_ph(pvars(PRES),pvars(ENTH),rhow_coeff,hw_coeff);
pvars(RHO_S) = rhos_ph(pvars(PRES),pvars(ENTH),rhos_coeff,hs_coeff);
pvars(H_W) = hw_p(pvars(PRES),hw_coeff);
pvars(H_S) = hs_p(pvars(PRES),hs_coeff);
pvars(SAT_S) = 1-pvars(SAT_W);
pvars(RHO) = pvars(RHO_W)*pvars(SAT_W) + pvars(SAT_S)*pvars(RHO_S);

   
% historical interpschem
if isfield(prop,'interpscheme')
    if strcmpi(prop.interpscheme,'actual') || strcmpi(prop.interpscheme,'fit')
        pvars(RHO) = interpRHO(pvars(PRES),pvars(ENTH),prop.rhodata,prop.interpscheme);
    else
end

end

