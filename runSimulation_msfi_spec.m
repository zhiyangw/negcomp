function [ converged,xk,xkHist ,Rnorm] = runSimulation_msfi_spec( xn,xk,prop,maxIter )

converged = 0;
iter = 1;
Rnorm = 0;
xkHist = [];

% 
% xf = xk;
% xk_hyp = zeros(prop.NB*2,1);
% for i = 1:prop.NB
%     xk_hyp(i) = rho_ph(xf(i),xf(i+prop.NB));
%     xk_hyp(i+prop.NB) = xf(i+prop.NB)*rho_ph(xf(i),xf(i+prop.NB));
% 
% end
% xkHist_hyp = xk_hyp;
outer_iter = 0;
while ~converged
    [massConverged,xk,mass_xkHist] = mass_loop(xn,xk,prop);
    xkHist = [xkHist mass_xkHist];
    if ~massConverged
        disp('mass failed')
        return 
    end
    xkHist = xk;
    [ converged,xk,xkHist_fim,Rnorm ] = runSimulation( xn,xk,prop,100,1);
    xkHist = [xkHist xkHist_fim];
    if converged
        break
    end
end

end

