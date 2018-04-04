function [prop,xn] = initialize(NB,dp,dh)


prop.dim = 0; % 0- dimensional 1-nondim
prop.NB = NB;
%prop.scenario = 'dh';
prop.chop =0;

prop.deps = 1e-4;
prop.perm = 1e-12;

length = 1;
prop.dx = length/prop.NB;
prop.dy = 1;
prop.dz = 1;

prop.V =prop.dx*prop.dy*prop.dz*ones(prop.NB,1);
prop.Reps = 1e-5;
prop.drho = 0.1;
prop.dE = 10;
prop.pinit = 8e6;
prop.Swinit = 0.1;
prop.hinit = h_pS(prop.pinit,prop.Swinit);
prop.poro = 0.15;
prop.rock = 0;
well.BHP = 9e6;
well.T = 1e-3;
well.blockNum = 1;
well.Hinj = 3.45e5;
well.type = 'bhp';
well2 = well;
well2.BHP = 1e6;
well2.T =0 ;
well2.blockNum = NB;
well2.type = 'bhp';
prop.wellList = {};
%prop.wellList(1) = well;
% no producer well
%prop.wellList(2) = well2;
prop.C_R = 4e6;
for i = 1:prop.NB-1
    prop.connList(i,:) = [i,i+1];
end
prop.primaryVariable = 'h';

xn1Init = pconvInv(prop.pinit,prop);
if strcmpi(prop.primaryVariable,'S')
    xn2Init = Sw(prop.pinit,prop.hinit);
elseif strcmpi(prop.primaryVariable,'h')
    xn2Init = hconvInv(prop.hinit,prop);
else
    error('unknown primary variable')
end

xn = [ones(prop.NB,1)*xn1Init;ones(prop.NB,1)*xn2Init];
%xn(3) = hconvInv(1e6,prop);
if exist('dp','var')
    prop.dp = dp;
else
    prop.dp = abs(xn(1)*prop.deps);
end
if exist('dh','var')
    prop.dh = dh;
else
    
    prop.dh = abs(xn(prop.NB+1)*prop.deps);
end
end


