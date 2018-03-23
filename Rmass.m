function [ R,Rnormed ] = Rmass( xn,xn1,prop ,block1Ind,block2Ind)

[PRES,ENTH,TEMP,SAT_W,SAT_S,RHO_W,RHO_S,RHO,H_W,H_S,VISC_W,VISC_S,KR_W,KR_S,NVAR] = offset();
% Calculate fluid properties for previous time step and current timestep
[pvars_1n] = fluidProperties(xn(block1Ind),xn(block1Ind+prop.NB),prop);
[pvars_2n] = fluidProperties(xn(block2Ind),xn(block2Ind+prop.NB),prop);
[pvars_1n1]= fluidProperties(xn1(block1Ind),xn1(block1Ind+prop.NB),prop);
[pvars_2n1] = fluidProperties(xn1(block2Ind),xn1(block2Ind+prop.NB),prop);

if pvars_1n1(PRES) > pvars_2n1(PRES)
    upstream_pvars = pvars_1n1;
else
    upstream_pvars = pvars_2n1;
end

R(1,1) = ((pvars_1n1(RHO)) ...
    - (pvars_1n(RHO))) ...
    + Tm(upstream_pvars,prop)*(prop.dt/prop.V(block1Ind))*(pvars_1n1(PRES) - pvars_2n1(PRES)); 

Rnormed(1,1)  = R(1,1) /(pvars_1n(RHO));
R(2,1) = (pvars_2n1(RHO)) ...
        - (pvars_2n(RHO)) ...
        + Tm(upstream_pvars,prop)*(prop.dt/prop.V(block2Ind))*(pvars_2n1(PRES) - pvars_1n1(PRES));


Rnormed(2,1)  = R(2,1) /(pvars_2n(RHO));
R = Rnormed;



end

