function [ J ] = J( xn,xn1,prop )
%Jacobian for FIM residual

J = zeros(prop.NB*2);
for i = 1:prop.NB
    dxn = xn1;
    dxn(i) = dxn(i) + prop.dp;
    J(:,i) = (R(xn,dxn,prop) - R(xn,xn1,prop))/prop.dp;
end
for i = 1:prop.NB
    dxn = xn1;
    dxn(i+prop.NB) = dxn(i+prop.NB) + prop.dh;
    J(:,i+prop.NB) = (R(xn,dxn,prop) - R(xn,xn1,prop))/(prop.dh);
end





end

