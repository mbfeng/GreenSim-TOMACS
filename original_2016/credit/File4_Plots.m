% Plots_Additional.m
% plots using estimators in SMC_ILR_MLR_Common and Additional_Estimators
close all

% CHANGE: load the data to be plotted
% load('results_10Kreps_justincase.mat')

%% figure 1: Current state plot
figure
hold on

errorbar(RMSE_Y_SMC, RMSE_E_SMC, '-', 'color',[0,204,204]/255, 'LineWidth',3)
errorbar(RMSE_Y_ILR, RMSE_E_ILR, '-.', 'color',[255,102,102]/255, 'LineWidth',3)
errorbar(RMSE_Y_WLR, RMSE_E_WLR, '-.', 'color',[255,204,0]/255, 'LineWidth',3)
errorbar(RMSE_Y_MLR, RMSE_E_MLR, '-.', 'color',[204,51,255]/255, 'LineWidth',3)
% errorbar(SKStart:NPeriods,RMSE_Y_SK1(SKStart:NPeriods),RMSE_E_SK1(SKStart:NPeriods),':', 'color',[0,0,204]/255, 'LineWidth',3)
% errorbar(RMSE_Y_Shr, RMSE_E_Shr, '-', 'color',[0.6,0.2,0], 'LineWidth',3)

set(gca,'FontSize',14,'FontWeight','bold','yscale','log','xscale','log','FontName','FixedWidth','xscale','log','yscale','log')
xlim([1,52])
ylim([1e-6,1e-3])
legend({'SMC','ILR','EWLR','MLR'}, 'FontSize',16, 'Location','southwest','FontWeight','bold')
legend boxoff
xlabel('Number of periods, n', 'FontSize',15,'FontWeight','bold')
title('period n variance for the current state','FontSize',15,'FontWeight','bold')
hold off

