

clear
close all
clc
pBound = [1e5 10e6];
x2Bound = [1e5 2.5e6];

tolerance = 0;
NB = 1;
nPlot = 10;
plotContours =1;
% intialize property values
[prop,xn] = initialize(NB,50,1);

p_hw = linspace(pBound(1),pBound(2),nPlot);
for i = 1:length(p_hw)
    hw_plot(i) = hw_p(p_hw(i));
end
dtList = [0.0045,0.12];
nPoints = 5;
nCont = 10;
for dtInd = 1:length(dtList)
    if dtInd == 1
        pBound = [2e5 1.2e6];
        x2Bound = [5e5 1e6];
    else
        pBound = [1e5 10e6];
        x2Bound = [1e5 2.5e6];
    end
    prop.dt =dtList(dtInd);
    pInit = linspace(pBound(1),pBound(2),nPoints);
    hInit = linspace(x2Bound(1),x2Bound(2),nPoints);
    [pCont,hCont] = meshgrid(linspace(pBound(1),pBound(2),nCont),linspace(x2Bound(1),x2Bound(2),nCont));
    
    figure(dtInd)
    
    [ converged,xk,xkHist,Rnorm ] = runSimulation( xn,xn,prop,100,0.05);
    plot(xkHist(1,:),xkHist(2,:),'k-','LineWidth',2.5);
    xlabel('Pressure (Pa)')
    ylabel('Enthalpy (J/kg)')
    hold on
    plot(xkHist(1,1),xkHist(2,1),'b.','markersize',20);
    plot(xkHist(1,end),xkHist(2,end),'r.','markersize',20);
    
    
    axis([pBound(1) pBound(2) x2Bound(1) x2Bound(2)]);
    disp('Plotting Contour')
    for pInd = 1:nCont
        disp(pInd)
        for hInd = 1:nCont
            xk = [pCont(pInd,hInd);hCont(pInd,hInd)];
            [Rout,~] = R(xn,xk,prop);
            Rcont(pInd,hInd) = log10(norm(Rout));
        end
    end
    contour(pCont,hCont,Rcont,20);
    colorbar
    disp('Plotting Newton Paths')
    for pInd = 1:nPoints
        disp(pInd)
        for hInd = 1:nPoints
            xInit = [pInit(pInd);hInit(hInd)];

            [ converged,xk,xkHist,Rnorm ] = runSimulation( xn,xInit,prop,350,0.05 );

            plot(xkHist(1,:),xkHist(2,:),'k-','LineWidth',1.5);
            xlabel('Pressure (Pa)')
            ylabel('Enthalpy (J/kg)')
            hold on
            plot(xkHist(1,1),xkHist(2,1),'k.','markersize',20);
            plot(xkHist(1,end),xkHist(2,end),'r.','markersize',20);
        end
    end
    figure(dtInd)
    plot(p_hw,hw_plot,'b-','linewidth',1.5);
    
    
end
