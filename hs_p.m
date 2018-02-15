function [ hs ] = hs_p( p ,h_coeff )

hconv = 1E4;
pconv = 10;
if exist('h_coeff','var')
    hs  = h_coeff(1)+h_coeff(2)./(p*pconv) ...
        +h_coeff(3).*(p*pconv).^(-2)...
        +h_coeff(4).*(p*pconv).^(2);
else
    hs  = 2.82282e10 - 3.91952e15./(p*pconv) ...
        + 2.54342e21.*(p*pconv).^(-2)...
        - 9.38879e-8.*(p*pconv).^(2);
end
hs  = hs /hconv;

end

