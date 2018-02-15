function [ p ] = pconv( ptilde,prop )

if prop.dim == 1
    p = ptilde*prop.BHP1 + prop.BHP1;
elseif prop.dim ==0
    p = ptilde;
end
end

