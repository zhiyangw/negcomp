

clear
close all
clc
pBound = [7.5e6 8.5e6];
x2Bound = [1e5 2.5e6];

tolerance = 0;
NB = 2;
nPlot = 10;
plotContours =1;
% intialize property values
[prop,xn] = initialize(NB,50,1);
p1 = 8e6;
p2 = 8.1e6;
xn = [p1;p2;h_pS(p1,0.9);h_pS(p2,0.5)];

prop.depth = [100,200];

p_hw = linspace(pBound(1),pBound(2),nPlot);
for i = 1:length(p_hw)
    hw_plot(i) = hw_p(p_hw(i));
end
dtList = [10000];
nPoints = 3;
nCont = 75;
prop.massMaxIter = 20;
prop.energyMaxIter = 20;
prop.damp = 1;
for dtInd = 1:length(dtList)
    % bounds to test the newton and contours
    prop.dt = dtList(dtInd);
    [ converged,xk,xkHist,Rnorm ] = runSimulation( xn,xn,prop,100,1);
    
end
figure

hold on
h(1) = plot(xkHist(1,:),xkHist(3,:),'k-','markersize',10,'linewidth',1);
plot(xkHist(1,1),xkHist(3,1),'r.','markersize',20,'linewidth',2);
plot(xkHist(1,end),xkHist(3,end),'b.','markersize',20,'linewidth',2);
h(2) = plot(xkHist(2,:),xkHist(4,:),'k-','markersize',10,'linewidth',1);
plot(xkHist(2,1),xkHist(4,1),'rd','markersize',10,'linewidth',2);
plot(xkHist(2,end),xkHist(4,end),'bd','markersize',10,'linewidth',2);
legend(h,'block1','block2')
plot(p_hw,hw_plot,'b-');
SwPlot = zeros(prop.NB,size(xkHist,2));
for i = 1:size(xkHist,2)
    SwPlot(1,i) = Sw(xkHist(1,i),xkHist(3,i));
    SwPlot(2,i) = Sw(xkHist(2,i),xkHist(4,i));
    
end
figure
hold on
plot(SwPlot(1,:),SwPlot(2,:),'k.-','markersize',10,'linewidth',1);
xlabel('block1')
ylabel('block2')
h(1) = plot(SwPlot(1,1),SwPlot(2,1),'r.','markersize',20,'linewidth',2);
h(2) = plot(SwPlot(1,end),SwPlot(2,end),'b.','markersize',20,'linewidth',2);


legend(h,'start','end')