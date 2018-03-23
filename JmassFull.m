function [ J ] = JmassFull( xn,xn1,prop )

J = zeros(prop.NB);
for i = 1:prop.NB
    dxn = xn1;
    dxn(i) = dxn(i) + prop.dp;
    J(:,i) = (RmassFull(xn,dxn,prop) - RmassFull(xn,xn1,prop))/prop.dp;
end


end

