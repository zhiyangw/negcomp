function [converged,xk,xkHist] = mass_loop(xn,xk,prop)
    pConv = 0;
    iter = 0;
    xkHist = xk;
    converged = 1;
    maxIter = prop.massMaxIter;
    while ~pConv
        [Rmass_out,Rnormed] = RmassFull(xn,xk,prop);
        
        if (iter ==0)
            disp('Mass')
            disp(norm(Rmass_out))
        end
        if norm(Rnormed) < prop.Reps
            if (iter == 0)
                pConv = 1;
            end
            break
        end
        dx = -JmassFull(xn,xk,prop)\Rmass_out;
        xk(1:length(dx)) = xk(1:length(dx)) + dx;
        %[prop,~,xk] =  updateX(xk,dx,prop.NB,prop.rhodata.coeff,0.1,prop,prop);
%         for i = 1:length(dx)
%             xk(i) = xk(i)+prop.damp*dx(i);
%         end
        %xk = updateSecondary_NB(xk_old,xk,prop);
        xkHist(:,size(xkHist,2) + 1) = xk;
        if any(xk < 0)
            disp('negative values')
            converged = 0;
            return
        end
        iter = iter + 1;
        if iter > maxIter
            disp('max iter')
            converged = 0;
            return
        end
        
    end
end