function [ p ] = pconvInv( ptilde,prop )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

if prop.dim == 1
    p = (ptilde-prop.BHP1)/prop.BHP1;
elseif prop.dim ==0
    p = ptilde;
end
end

