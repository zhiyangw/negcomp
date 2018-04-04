function [Te] = Tes(pvars,prop)
[PRES,ENTH,TEMP,SAT_W,SAT_S,RHO_W,RHO_S,RHO,H_W,H_S,VISC_W,VISC_S,KR_W,KR_S,NVAR] = offset();
Te = (prop.perm*prop.dy*prop.dz/prop.dx)*...
    (pvars(H_S)*pvars(RHO_S)*pvars(KR_S)/pvars(VISC_S));

end