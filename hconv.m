function [ h ] = hconv( htilde,prop )

if prop.dim == 1
    h = htilde*prop.hinit;
elseif prop.dim ==0
    h = htilde;
end


end

