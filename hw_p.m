function [ hw ] = hw_p( p ,h_coeff )

hconv = 1E4;
pconv = 10;
if exist('h_coeff','var')
    hw = h_coeff(1) + h_coeff(2)*(p*pconv) +h_coeff(3)*(p*pconv).^2 	+ ...
        h_coeff(4)*(p*pconv).^3			+h_coeff(5)./(p*pconv) ...
        +h_coeff(6).*(p*pconv).^(-2)		+h_coeff(7).*(p*pconv).^(-3);
else
    hw = 7.30984e9 + 1.29239e2*(p*pconv) -1.00333e-6*(p*pconv).^2 	+ ...
        3.9881e-15*(p*pconv).^3			- 9.09697e15./(p*pconv) ...
        +1.29267e22.*(p*pconv).^(-2)		-6.28359e27.*(p*pconv).^(-3);
end

hw  = hw /hconv;

end

