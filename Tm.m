function [Tm] = Tm(pvars,prop)
[PRES,ENTH,TEMP,SAT_W,SAT_S,RHO_W,RHO_S,RHO,H_W,H_S,VISC_W,VISC_S,KR_W,KR_S,W,DWDP,DWDN,DWDV,NVAR]  = offset();
Tm = (prop.perm*prop.dy*prop.dz/prop.dx)*...
    (pvars(RHO_W)*pvars(KR_W)/pvars(VISC_W) + pvars(RHO_S)*pvars(KR_S)/pvars(VISC_S));

end