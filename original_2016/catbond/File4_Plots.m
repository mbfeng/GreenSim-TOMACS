% Plots_Additional.m
% plots using estimators in SMC_ILR_MLR_Common and Additional_Estimators
close all

% CHANGE: load the data to be plotted
% load('results_10Kreps_justincase.mat')

%% figure 1: Current state plot
figure
hold on

ref_x = logspace(0,log(100)/log(10),20);
ref_y =(ref_x.^(-1)*RMSE_Y_SMC(1));

plot(1:NPeriods, RMSE_Y_SMC, '-', 'color',[0,204,204]/255, 'LineWidth',4)
plot(1:NPeriods, RMSE_Y_ILR, '-.', 'color',[255,102,102]/255, 'LineWidth',4)
plot(1:NPeriods, RMSE_Y_WLR, '-.', 'color',[255,204,0]/255, 'LineWidth',4)
plot(1:NPeriods, RMSE_Y_MLR, '-.', 'color',[204,51,255]/255, 'LineWidth',4)
% plot(SKStart:NPeriods,RMSE_Y_SK1(SKStart:NPeriods),':', 'color',[0,0,204]/255, 'LineWidth',4)
plot(ref_x,ref_y, '-', 'LineWidth',4, 'color',[0,0,0]/255)
% plot(2:NPeriods, RMSE_Y_Shr(2:NPeriods), '-', 'color',[0.6,0.2,0], 'LineWidth',4)

% xlim([0,100])
% ylim([1e-2,2e1])

legend({'SMC','ILR','EWLR','MLR','Slope -1'}, 'FontSize',16, 'Location','southwest','FontWeight','bold')
legend boxoff
set(gca,'FontSize',14,'FontWeight','bold','yscale','log','xscale','log','FontName','FixedWidth','xscale','log','yscale','log')
xlabel('Number of periods, n', 'FontSize',15,'FontWeight','bold')
title('period n MSE for the current state','Position',[10 1e-3],'FontSize',15,'FontWeight','bold')
hold off

%% figure 2: Fixed state plot (just mi- and hi-states in paper)
figure
ha = tight_subplot(1,2,[.1 .05],[.1 .05],[.05 .05]);

axes(ha(1));
ref_x = logspace(0,log(100)/log(10),20);
ref_y = (ref_x.^(-1)*RMSE_Y_SMC_mi(1));

hold on
plot(RMSE_Y_SMC_mi, '-',  'color',[0,204,204]/255, 'LineWidth',4)
plot(RMSE_Y_ILR_mi, '-.', 'color',[255,102,102]/255, 'LineWidth',4)
plot(RMSE_Y_WLR_mi, '-.', 'color',[255,204,0]/255, 'LineWidth',4)
plot(RMSE_Y_MLR_mi, '-.', 'color',[204,51,255]/255, 'LineWidth',4)
% plot(SKStart:NPeriods,RMSE_Y_SK1_mi(SKStart:end), ':', 'color',[0,0,204]/255, 'LineWidth',4)
plot(ref_x,ref_y, '-', 'LineWidth',4, 'color',[0,0,0]/255)
% plot(2:NPeriods, RMSE_Y_Shr_mi(2:NPeriods), '-', 'color',[0.6,0.2,0], 'LineWidth',4)

set(gca,'FontSize',14,'FontWeight','bold','FontName','FixedWidth', ...
    'yscale','log','xscale','log', ...
    'YTickLabelMode','auto','XTickLabelMode','auto')
legend({'SMC','ILR','EWLR','MLR','Slope -1'}, 'FontSize',16, 'Location','southwest','FontWeight','bold')
legend boxoff
xlabel('Number of periods, n', 'FontSize',15,'FontWeight','bold')
title('period n MSE for center state-x_{mi}','Position',[10 2e-2], 'FontSize',15,'FontWeight','bold')
xlim([0,100])
ylim([3e-6,2e-2])
hold off

axes(ha(2));
ref_x = logspace(0,log(100)/log(10),20);
ref_y = (ref_x.^(-1)*RMSE_Y_SMC_hi(1));

hold on
plot(RMSE_Y_SMC_hi, '-',  'color',[0,204,204]/255, 'LineWidth',4)
plot(RMSE_Y_ILR_hi, '-.', 'color',[255,102,102]/255, 'LineWidth',4)
plot(RMSE_Y_WLR_hi, '-.', 'color',[255,204,0]/255, 'LineWidth',4)
plot(RMSE_Y_MLR_hi, '-.', 'color',[204,51,255]/255, 'LineWidth',4)
% plot(SKStart:NPeriods, RMSE_Y_SK1_hi(SKStart:end), ':', 'color',[0,0,204]/255, 'LineWidth',4)
plot(ref_x,ref_y, '-', 'LineWidth',4, 'color',[0,0,0]/255)
% plot(2:NPeriods, RMSE_Y_Shr_hi(2:NPeriods), '-', 'color',[0.6,0.2,0], 'LineWidth',4)

set(gca,'FontSize',14,'FontWeight','bold','FontName','FixedWidth', ...
    'yscale','log','xscale','log', ...
    'YTickLabelMode','auto','XTickLabelMode','auto')
legend({'SMC','ILR','EWLR','MLR','Slope -1'}, 'FontSize',16, 'Location','southwest','FontWeight','bold')
legend boxoff
xlabel('Number of periods, n', 'FontSize',15,'FontWeight','bold')
title('period n MSE for extreme state-x_{hi}','Position',[10 2e-2],'FontSize',15,'FontWeight','bold')
xlim([0,100])
ylim([3e-6,2e-2])
hold off

% % %% figure 3: Histogram of Parameters
% % figure
% % ha = tight_subplot(1,2,[.1 .05],[.1 .05],[.05 .05]);
% % 
% % axes(ha(1));
% % histogram(Paths(:,1),200,'EdgeAlpha',0,'Normalization','pdf')
% % 
% % axes(ha(2));
% % histogram(Paths(:,2),200,'EdgeAlpha',0,'Normalization','pdf')