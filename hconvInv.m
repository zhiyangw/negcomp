function [ h ] = hconvInv( htilde,prop )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
h = htilde/prop.hinit;
if prop.dim == 1
    h = htilde/prop.hinit;
elseif prop.dim ==0
    h = htilde;
end
end

