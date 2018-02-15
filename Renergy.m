function [ R,Rnormed ] = Renergy( xn,xnf,xn1,prop ,block1Ind,block2Ind)

[PRES,ENTH,TEMP,SAT_W,SAT_S,RHO_W,RHO_S,RHO,H_W,H_S,VISC_W,VISC_S,KR_W,KR_S,NVAR] = offset();


% xn = [p1 p2 h1 h2]
[pvars_1n] = fluidProperties(xn(block1Ind),xn(block1Ind+prop.NB),prop);
[pvars_2n] = fluidProperties(xn(block2Ind),xn(block2Ind+prop.NB),prop);


[pvars_1n1f]= fluidProperties(xnf(block1Ind),xnf(block1Ind+prop.NB),prop);
[pvars_2n1f] = fluidProperties(xnf(block2Ind),xnf(block2Ind+prop.NB),prop);

%[pvars_1n1]= fluidProperties_rhoh(pvars_1n1f(RHO),xn1(i+prop.NB),prop);
%[pvars_2n1] = fluidProperties_rhoh(pvars_2n1f(RHO),xn1(j+prop.NB),prop);
[pvars_1n1]= fluidProperties(xn1(block1Ind),xn1(block1Ind+prop.NB),prop);
[pvars_2n1] = fluidProperties(xn1(block2Ind),xn1(block2Ind+prop.NB),prop);

if pvars_1n1(PRES) > pvars_2n1(PRES)
    upstream_pvars = pvars_1n1f;
else
    upstream_pvars = pvars_2n1f;
end

R(1,1) = (pvars_1n1f(RHO)*pvars_1n1(ENTH) ...
    - pvars_1n(RHO)*pvars_1n(ENTH) + (0.8/0.2)*prop.C_R*(pvars_1n1(TEMP)-pvars_1n(TEMP))) ...
    + Te(upstream_pvars,prop)*(prop.dt/prop.V(block1Ind))*(pvars_1n1f(PRES) - pvars_2n1f(PRES));...



Rnormed(1,1)  = R(1,1) /((pvars_1n(RHO)*pvars_1n(ENTH)));

R(2,1) = (pvars_2n1f(RHO)*pvars_2n1(ENTH) ...
    - pvars_2n(RHO)*pvars_2n(ENTH) + (0.8/0.2)*0*(pvars_2n1(TEMP)-pvars_2n(TEMP))) ...
    + Te(upstream_pvars,prop)*(prop.dt/prop.V(block2Ind))*(pvars_2n1f(PRES) - pvars_1n1f(PRES));...
Rnormed(2,1)  = R(2,1) /((pvars_2n(RHO))*pvars_2n(ENTH));
    
R = Rnormed;    
    

end

