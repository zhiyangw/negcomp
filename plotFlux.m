
clear all; close all
nS = 100;
h_plot = linspace(5e5,2.9e6,nS);
S_plot = h_plot;%linspace(0,1,nS);
nP = 3;
dS = S_plot(2)-S_plot(1);
pList = linspace(1,10,nP)*1e6;
Ng_list = [-20,0,20];
figure('pos',[10 10 900 1200])
figure('pos',[10 10 900 1200])
figure('pos',[10 10 900 1200])
figure('pos',[10 10 900 1200])

for pInd = 1:length(pList)
    pLegend = {};
    p = pList(pInd);
    for NgInd = 1:length(Ng_list)
        Ng = Ng_list(NgInd);
        coeff = [1.00207; 4.42607; -5.47456; 5.02875; -1.24791; ...
            -2.26162; 4.38441; -1.79088; 3.69276; 5.17644;...
            7.30984; 1.29239; -1.00333; 3.9881; -9.09697; 1.29267; -6.28359; ...
            2.82282; -3.91952; 2.54342; -9.38879];
        for i = 1:nS
            h_val = h_plot(i); %h_pS(p,S_plot(i));
            Sval = Sw(p,h_val); %S_plot(i);

            fw(i) = (krw(p,h_val,coeff)/viscw(p,h_val,coeff))/(krw(p,h_val,coeff)/viscw(p,h_val,coeff)+krs(p,h_val,coeff)/viscs(p,h_val,coeff));
            fw(i) = fw(i) - krs(p,h_val,coeff)*fw(i)*Ng;
            Fm(i) = rhow_p(p)*Sval*fw(i)+rhos_p(p)*(1-Sval)*(1-fw(i));
            Fe(i) = rhow_p(p)*hw_p(p)*Sval*fw(i)+rhos_p(p)*hs_p(p)*(1-Sval)*(1-fw(i));
        end
        Fm = Fm/(rhow_p(p)-rhos_p(p));
        Fe = Fe/(rhow_p(p)*hw_p(p)-rhos_p(p)*hs_p(p));
        figure(1)
        subplot(3,2,(pInd-1)*2+1)
        plot(S_plot,Fm)
        hold on
        subplot(3,2,(pInd-1)*2+2)
        plot(S_plot,Fe)
        hold on
        figure(2)
        mid_sw =0.5*(S_plot(1:end-1)+S_plot(2:end));
        dFm = diff(Fm)/dS;
        dFe = diff(Fe)/dS;
        subplot(3,2,(pInd-1)*2+1)
        plot(mid_sw,dFm)
        hold on
        subplot(3,2,(pInd-1)*2+2)
        plot(mid_sw,dFe)
        hold on
        figure(3)
        dmidsw = mid_sw(2) - mid_sw(1);
        midmid_sw =0.5*(mid_sw(1:end-1)+mid_sw(2:end));
        ddFm = diff(dFm)/dmidsw;
        ddFe = diff(dFe)/dmidsw;
        subplot(3,2,(pInd-1)*2+1)
        plot(midmid_sw,ddFm)
        hold on
        subplot(3,2,(pInd-1)*2+2)
        plot(midmid_sw,ddFe)
        hold on
        figure(4)
        subplot(3,2,(pInd-1)*2+1)
        plot(S_plot,fw)
        hold on
        subplot(3,2,(pInd-1)*2+2)
        plot(mid_sw,diff(fw)/dS)
        hold on
        pLegend{end+1} = ['Ng = ' num2str(Ng)]
    end
    figure(1)
    subplot(3,2,(pInd-1)*2+1)
    plot([min(S_plot) max(S_plot)],[0 0],'k--')
    title(['p = ' num2str(p)])
    xlabel('Sw')
    ylabel('Fm')
    legend(pLegend,'location','northwest')
    subplot(3,2,(pInd-1)*2+2)
    plot([min(S_plot) max(S_plot)],[0 0],'k--')
    title(['p = ' num2str(p)])
    xlabel('Sw')
    ylabel('Fe')
    legend(pLegend,'location','northwest')
    figure(2)
    subplot(3,2,(pInd-1)*2+1)
    plot([min(S_plot) max(S_plot)],[0 0],'k--')
    title(['p = ' num2str(p)])
    xlabel('Sw')
    ylabel('dFm')
    legend(pLegend,'location','northwest')
    subplot(3,2,(pInd-1)*2+2)
    plot([min(S_plot) max(S_plot)],[0 0],'k--')
    title(['p = ' num2str(p)])
    xlabel('Sw')
    ylabel('dFe')
    legend(pLegend,'location','northwest')
    figure(3)
    subplot(3,2,(pInd-1)*2+1)
    plot([min(S_plot) max(S_plot)],[0 0],'k--')
    title(['p = ' num2str(p)])
    xlabel('Sw')
    ylabel('d2Fm')
    legend(pLegend,'location','northwest')
    subplot(3,2,(pInd-1)*2+2)
    plot([min(S_plot) max(S_plot)],[0 0],'k--')
    title(['p = ' num2str(p)])
    xlabel('Sw')
    ylabel('d2Fe')
    figure(4)
    subplot(3,2,(pInd-1)*2+1)
    plot([min(S_plot) max(S_plot)],[0 0],'k--')
    title(['p = ' num2str(p)])
    xlabel('Sw')
    ylabel('fw')
    legend(pLegend,'location','northwest')
    subplot(3,2,(pInd-1)*2+2)
    plot([min(S_plot) max(S_plot)],[0 0],'k--')
    title(['p = ' num2str(p)])
    xlabel('Sw')
    ylabel('dfw')
    legend(pLegend,'location','northwest')
end

figure(1)
saveToPdf('flux_plot')
figure(2)
saveToPdf('dfds_plot')
figure(3)
saveToPdf('d2fds2_plot')
figure(4)
saveToPdf('fw_plot')
