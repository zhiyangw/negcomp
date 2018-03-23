function [ J ] = JenergyFull( xn,xk1,xn1,prop )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
J = zeros(prop.NB);
for i = 1:prop.NB
    dxn = xn1;
    dxn(i+prop.NB) = dxn(i+prop.NB) + prop.dh;
    J(:,i) = (RenergyFull(xn,xk1,dxn,prop) - RenergyFull(xn,xk1,xn1,prop))/prop.dh;
end

end

