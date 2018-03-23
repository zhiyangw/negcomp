function [ converged, xk, xkHist ] = energy_loop( xn,xk,prop )

    hConv = 0;
    iter = 0;
    converged = 1;
    xkHist = xk;
    maxIter = prop.energyMaxIter;
    while ~hConv
        Rh_out = RenergyFull(xn,xk,xk,prop);
        
        if (iter ==0)
            disp('Energy')
            disp(norm(Rh_out))
        end
        if norm(Rh_out) < prop.Reps
            if (iter == 0)
                hConv = 1;
            end
            break
        end
        dx = -JenergyFull( xn,xk,xk,prop )\Rh_out;
        for i = 1:length(dx)
            xk(i+prop.NB) = xk(i+prop.NB)+prop.damp*dx(i);
        end
        %xk = updateSecondary_NB(xk_old,xk,prop);
        xkHist(:,size(xkHist,2) + 1) = xk;
        if any(xk < 0) 
            disp('negative values')
            converged = 0;
            return
        end
        if any(xk > 2e7) 
            disp('too large enthalpy')
            converged = 0;
            return
        end
%         for i = 1:prop.NB
%             if getStatus(xk(i:prop.NB:end),prop) == 0
%                 %                 pConv = 1;
%                 %                 break
%             end
%         end
        iter = iter + 1;
        if iter > maxIter
            break
        end
        
    end

end

