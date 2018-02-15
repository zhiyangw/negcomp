function [ converged,xk,xkHist ,Rnorm] = runSimulation( xn,xk,prop,maxIter,chop )

    converged = 0;
    xkHist = [xk];
    iter = 1;
    
    while ~converged
        
        [Rout,~] = R(xn,xk,prop);
        dx = -J(xn,xk,prop)\Rout;
        xkHist(:,iter) = xk;
                
        xk = xk+chop*dx;
        
        if (norm(Rout) < prop.Reps)
            converged = 1;
        end
        if iter > maxIter
            break
        end
        if any(xk < 0)
            xkHist(:,iter+1) = xk;
            break
        end
        iter = iter+1;
    end
    Rnorm = norm(Rout);
end

