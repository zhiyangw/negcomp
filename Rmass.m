function [ R,Rnormed ] = Rmass( xn,xn1,prop ,block1Ind,block2Ind)

[PRES,ENTH,TEMP,SAT_W,SAT_S,RHO_W,RHO_S,RHO,H_W,H_S,VISC_W,VISC_S,KR_W,KR_S,NVAR] = offset();
% Calculate fluid properties for previous time step and current timestep
[pvars_1n] = fluidProperties(xn(block1Ind),xn(block1Ind+prop.NB),prop);
[pvars_2n] = fluidProperties(xn(block2Ind),xn(block2Ind+prop.NB),prop);
[pvars_1n1]= fluidProperties(xn1(block1Ind),xn1(block1Ind+prop.NB),prop);
[pvars_2n1] = fluidProperties(xn1(block2Ind),xn1(block2Ind+prop.NB),prop);

if pvars_1n1(PRES) > pvars_2n1(PRES)
    upw_pvars = pvars_1n1;
    ups_pvars = pvars_1n1;
else
    upw_pvars = pvars_2n1;
    ups_pvars = pvars_2n1;
end

R(1,1) = ((pvars_1n1(RHO)) ...
    - (pvars_1n(RHO)));

R(2,1) = (pvars_2n1(RHO)) ...
    - (pvars_2n(RHO));

potdifw = (pvars_1n1(PRES)-(prop.depth(block1Ind)*pvars_1n1(RHO_W)))...
    - (pvars_2n1(PRES)-(prop.depth(block2Ind)*pvars_2n1(RHO_W)));
if potdifw >0
    upw_pvars = pvars_1n1;
else
    upw_pvars = pvars_2n1;
end
R(1,1) = R(1,1) + (prop.dt/prop.V(block1Ind))*(Tmw(upw_pvars,prop))*potdifw;
R(2,1) = R(2,1) - (prop.dt/prop.V(block1Ind))*(Tmw(upw_pvars,prop))*potdifw;


potdifs = (pvars_1n1(PRES)-(prop.depth(block1Ind)*9.81*pvars_1n1(RHO_S)))...
    - (pvars_2n1(PRES)-(prop.depth(block2Ind)*9.81*pvars_2n1(RHO_S)));
if potdifs >0
    ups_pvars = pvars_1n1;
else
    ups_pvars = pvars_2n1;
end
R(1,1) = R(1,1) + (prop.dt/prop.V(block1Ind))*(Tms(ups_pvars,prop))*potdifs;
R(2,1) = R(2,1) - (prop.dt/prop.V(block1Ind))*(Tms(ups_pvars,prop))*potdifs;


Rnormed(1,1)  = R(1,1) /(pvars_1n(RHO));
Rnormed(2,1)  = R(2,1) /(pvars_2n(RHO));

R = Rnormed;



end

