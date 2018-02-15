function [ h ] = h_pS( p,S_w,coeff )

%rho_w(5), rho_s(5), h_w(7), h_s(4)
if exist('coeff','var')
    [ rhow_coeff,rhos_coeff,hw_coeff ,hs_coeff] = coeff_to_fun_coeff( coeff );
    % Calculate actual Sw

    rhow= rhow_p(p,rhow_coeff);
    rhos= rhos_p(p,rhos_coeff);

    hw= hw_p(p,hw_coeff);
    hs= hs_p(p,hs_coeff);
    S_s = 1-S_w;
else
    rhow= rhow_p(p);
    rhos= rhos_p(p);
    hw= hw_p(p);
    hs= hs_p(p);
    S_s = 1-S_w;
    
end
h = (hw.*rhow.*S_w + hs.*rhos.*S_s)./(rhow.*S_w + rhos.*S_s);


end

