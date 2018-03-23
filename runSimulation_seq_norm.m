function [ converged,xk,xkHist ,Rnorm] = runSimulation_seq_norm( xn,xk,prop,maxIter )

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
    
    
    [ energyConverged, xk, energy_xkHist ] = energy_loop( xn,xk,prop );
    xkHist = [xkHist energy_xkHist];
    if ~energyConverged
        disp('energy failed')
        return 
    end 
    
    Rmass_out = RmassFull(xn,xk,prop);
    Rh_out = RenergyFull(xn,xk,xk,prop);
    if (norm(Rmass_out) <prop.Reps && norm(Rh_out) < prop.Reps && energyConverged )
        converged = 1;
        break
    end
    outer_iter = outer_iter +1;
    disp(outer_iter);
    if outer_iter > maxIter
        disp('max outer loops')
        return
    end
end

end

