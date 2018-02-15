function [ R,Rnormed ] = RmassFull( xn,xn1,prop )

[PRES,ENTH,TEMP,SAT_W,SAT_S,RHO_W,RHO_S,RHO,H_W,H_S,VISC_W,VISC_S,KR_W,KR_S,NVAR] = offset();
R = zeros(prop.NB,1);
Rnormed = zeros(prop.NB,1);
% block transmissibilities
if prop.NB > 1
    % This probably currently doesn't work have to fix with the new shiz
    for i = 1:size(prop.connList,1)
        blocki = prop.connList(i,1);
        blockj = prop.connList(i,2);
        [RBlock,RnormedBlock] = Rmass(xn,xn1,prop,blocki,blockj);
        R([blocki,blockj]) = R([blocki,blockj]) + RBlock;
        Rnormed([blocki,blockj]) = Rnormed([blocki,blockj]) + RnormedBlock;
    end
else
    [RBlock,RnormedBlock] = Rmass(xn,xn1,prop,1,1);
    R = RBlock(1);
    Rnormed = RnormedBlock(1);
end
% add wells
for i = 1:size(prop.wellList,2)
    well = prop.wellList(i);
    if strcmp(well.type,'bhp')
        [pvars] = fluidProperties(xn(well.blockNum),xn(well.blockNum+prop.NB),prop);
        [pvarsn1] = fluidProperties(xn1(well.blockNum),xn1(well.blockNum+prop.NB),prop);
        R(well.blockNum) = R(well.blockNum) - well.T*(prop.dt/prop.V(well.blockNum))*(well.BHP-pvarsn1(PRES))...
            /(pvars(RHO));
        Rnormed(well.blockNum) = Rnormed(well.blockNum) - well.T*(prop.dt/prop.V(well.blockNum))*(well.BHP-pvarsn1(PRES))...
            /(pvars(RHO));
    elseif strcmp(well.type,'rate')
        [pvars] = fluidProperties(xn(well.blockNum),xn(well.blockNum+prop.NB),prop);
        R(well.blockNum) = R(well.blockNum) - well.rate*(prop.dt/prop.V)/(pvars(RHO));
    end
end


end

