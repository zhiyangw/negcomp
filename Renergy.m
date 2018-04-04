function [ R,Rnormed ] = Renergy( xn,xnf,xn1,prop ,block1Ind,block2Ind)

[PRES,ENTH,TEMP,SAT_W,SAT_S,RHO_W,RHO_S,RHO,H_W,H_S,VISC_W,VISC_S,KR_W,KR_S,NVAR] = offset();


% xn = [p1 p2 h1 h2]
[pvars_1n] = fluidProperties(xn(block1Ind),xn(block1Ind+prop.NB),prop);
[pvars_2n] = fluidProperties(xn(block2Ind),xn(block2Ind+prop.NB),prop);


[pvars_1n1f]= fluidProperties(xnf(block1Ind),xnf(block1Ind+prop.NB),prop);
[pvars_2n1f] = fluidProperties(xnf(block2Ind),xnf(block2Ind+prop.NB),prop);

[pvars_1n1]= fluidProperties(xn1(block1Ind),xn1(block1Ind+prop.NB),prop);
[pvars_2n1] = fluidProperties(xn1(block2Ind),xn1(block2Ind+prop.NB),prop);

if pvars_1n1(PRES) > pvars_2n1(PRES)
    upw_pvars = pvars_1n1;
    ups_pvars = pvars_1n1;
else
    upw_pvars = pvars_2n1;
    ups_pvars = pvars_2n1;
end

R(1,1) = (pvars_1n1f(RHO)*pvars_1n1(ENTH) ...
    - pvars_1n(RHO)*pvars_1n(ENTH) + (0.8/0.2)*prop.C_R*(pvars_1n1(TEMP)-pvars_1n(TEMP)));

R(2,1) = (pvars_2n1f(RHO)*pvars_2n1(ENTH) ...
    - pvars_2n(RHO)*pvars_2n(ENTH) + (0.8/0.2)*prop.C_R*(pvars_2n1(TEMP)-pvars_2n(TEMP)));

  
potdifw = (pvars_1n1(PRES)-(prop.depth(block1Ind)*pvars_1n1(RHO_W)))...
    - (pvars_2n1(PRES)-(prop.depth(block2Ind)*pvars_2n1(RHO_W)));
if potdifw >0
    upw_pvars = pvars_1n1;
else
    upw_pvars = pvars_2n1;
end

R(1,1) = R(1,1) + (prop.dt/prop.V(block1Ind))*(Tew(upw_pvars,prop))*potdifw;
R(2,1) = R(2,1) - (prop.dt/prop.V(block2Ind))*(Tew(upw_pvars,prop))*potdifw;

potdifs = (pvars_1n1(PRES)-(prop.depth(block1Ind)*9.81*pvars_1n1(RHO_S)))...
    - (pvars_2n1(PRES)-(prop.depth(block2Ind)*9.81*pvars_2n1(RHO_S)));
if potdifs >0
    ups_pvars = pvars_1n1;
else
    ups_pvars = pvars_2n1;
end
R(1,1) = R(1,1) + (prop.dt/prop.V(block1Ind))*(Tes(ups_pvars,prop))*potdifs;
R(2,1) = R(2,1) - (prop.dt/prop.V(block2Ind))*(Tes(ups_pvars,prop))*potdifs;


Rnormed(1,1)  = R(1,1) /((pvars_1n(RHO)*pvars_1n(ENTH)));    
Rnormed(2,1)  = R(2,1) /((pvars_2n(RHO)*pvars_2n(ENTH)));

R = Rnormed;    
    

end

