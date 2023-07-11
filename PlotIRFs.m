function PlotIRFs(Results,figureName)
%PlotIRFs Replication of figures with impulse responses in paper.
%   Ewoud Quaghebeur, July 2023

T=size(Results.REE_IRF_eg.y,1);

variablesMAT={
    'y'     'Output ($Y$)'
    'c'     'Consumption ($C$)'
    'i'     'Investment ($I$)'
    'n'     'Hours worked ($N$)'
    'w'     'Real wages ($W$)'
    'r'     'Nominal interest rate ($R$)'
    'rk'    'Rental rate of capital ($r^k$)'
    };

figure('Name',figureName);
for iPlot=1:7
    subplot(4,2,iPlot);
    plot(1:T,100*Results.REE_IRF_eg.y(:,strcmp(Results.REE_IRF_eg.names,variablesMAT(iPlot,1))),...
        1:T,100*Results.RL_eg.y(strcmp(Results.REE_IRF_eg.names,variablesMAT(iPlot,1)),:),'--',...
        1:T,zeros(T),'Color','k');
    title(variablesMAT(iPlot,2),'interpreter','latex');
end
subplot(4,2,8);
plot(1:T,100*Results.REE_IRF_eg.y(:,strcmp(Results.REE_IRF_eg.names,'g')),...
    1:T,100*Results.RL_eg.w(strcmp(Results.RL_eg.names_z,'g'),:),'--',...
    1:T,zeros(T),'Color','k');
title('Government spending ($G$)','interpreter','latex');
legend('Rational Expectations','Adaptive Learning')
end