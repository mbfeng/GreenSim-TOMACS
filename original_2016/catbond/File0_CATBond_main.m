% File0_CATBond_main.m
% Green simulation for CAT bond example
% Oct 1, 2016
% -- Four green estimators with common random numbers (CRN)
%    1. Individual Likelihood Ratio (ILR)
%    2. Weighted Likelihood Ratio (WLR)
%    3. Multiple Likelihood Ratio (MLR)
%    4, Green Stochastic Kriging (GSK)

%% Workspace initialization
clear all   % clear all variable in workspace
close all   % close all figures in current session
rng(1);     % fix random seed for reproducibility
clc;        % clear current screen
addpath([pwd, '\StochasticKriging']) % load necessary functions for GSK

%% Experiment Parameters
% ----- Experiment design parameters
NPaths = 1e4;   % number of sample paths (1e2 for speed test, 1e4 in paper)
NPeriods = 1e2;   % number of periods per sample path
NReps = 1e2;    % number of replications per period
SKStart = 5;    % starting period of stochastic kriging

% ----- CAT bond parameters
K = 16;     % triggering point of Cat bond
p = 0.5;    % recovery rate in Cat bond

% ----- State process parameters
% xk+1 = mu + eta * xk + eps, var[eps] = sigma^2
% AR(1) process <===> Euler discretization of OU-process
mu = [0.0, 0.0]; % constant long-run mean
eta = [0.6, 0.5]; % autoregressive parameter
sigma = [0.8, 0.5]; % standard deviation of noise term
station_mean = mu./(1 - eta);   % stationary mean
station_std = sqrt((sigma.^2)./(1 - (eta).^2));   % stationary std

% ----- State transformation
% bounds for Poisson mean (lambda) & exponential mean (theta)
lambda_ub = 4; lambda_lb = 2;
theta_ub = 6; theta_lb = 4;
% sigmoid function with specified bounds
sigmoid = @(lb, ub, x) lb + (ub-lb)./(1 + exp(x));

% ----- Corner and certer parameters in the numerical experiment
lambda_lo = lambda_lb; theta_lo = theta_lb;    % bottom-left corner of parameter space
lambda_hi = lambda_ub; theta_hi = theta_ub;    % top-right corner of parameter space

% center of parameter space
lambda_mi = 0.5 * (lambda_lo + lambda_hi);
theta_mi = 0.5 * (theta_lo + theta_hi);

%% Simulate sample paths
% ----- placeholders for sample paths and true prices
Paths = zeros(NPeriods*NPaths, 2);
PriceTrue = zeros(NPeriods*NPaths, 1);

path = 0;
for i = 1:NPeriods*NPaths
    if(mod(i, NPeriods) == 1)
        path = path + 1
        % every paths starts from stationary distribution
        Paths(i,:) = normrnd(station_mean, station_std);
    else
        % sample paths evolves as an AR(1) process
        Paths(i,:) = mu + eta.*Paths(i-1,:) + normrnd([0, 0], sigma);
    end
end

% transformation from state space to parameter space
Paths(:,1) = sigmoid(lambda_lb, lambda_ub, Paths(:,1));
Paths(:,2) = sigmoid(theta_lb, theta_ub, Paths(:,2));

% ----- Calculate true prices in each period and each sample path
NTerms = 40; % number of approximation terms in calculating true prices
for i = 1:NPeriods*NPaths
    PriceTrue(i) = CatBondPrice(Paths(i,1), Paths(i,2), K, p, NTerms);
end

Price_lo = CatBondPrice(lambda_lo, theta_lo, K, p, NTerms);
Price_mi = CatBondPrice(lambda_mi, theta_mi, K, p, NTerms);
Price_hi = CatBondPrice(lambda_hi, theta_hi, K, p, NTerms);

%% Simulation model for CAT Bond pricing
% ----- placeholders for results
PriceSMC = zeros(NPeriods*NPaths, 1);   % Standard Monte Carlo (SMC)
PriceILR = zeros(NPeriods*NPaths, 1);   % Individual Likelihood Ratio (ILR)
PriceWLR = zeros(NPeriods*NPaths, 1);   % Weighted Likelihood Ratio (WLR)
PriceMLR = zeros(NPeriods*NPaths, 1);   % Mixture Likelihood Ratio (MLR)
PriceGSK = zeros(NPeriods*NPaths, 1);   % Green Stochastic Kriging (GSK)

PriceILR_mi = zeros(NPeriods*NPaths, 1);    % ILR estimate for mi-params
PriceILR_hi = zeros(NPeriods*NPaths, 1);    % ILR estimate for hi-params
PriceWLR_mi = zeros(NPeriods*NPaths, 1);    % WLR estimate for mi-params
PriceWLR_hi = zeros(NPeriods*NPaths, 1);    % WLR estimate for hi-params
PriceMLR_mi = zeros(NPeriods*NPaths, 1);    % MLR estimate for mi-params
PriceMLR_hi = zeros(NPeriods*NPaths, 1);    % MLR estimate for hi-params
PriceGSK_mi = zeros(NPeriods*NPaths, 1);    % GSK estimate for mi-params
PriceGSK_hi = zeros(NPeriods*NPaths, 1);    % GSK estimate for hi-params

% ----- main simulation model
for path = 1:NPaths
    tic
    NT = zeros(NPeriods*NReps, 1);  % PILRson r.v.
    CT = zeros(NPeriods*NReps, 1);  % Gamma r.v.
    IT = zeros(NPeriods*NReps, 1);  % Indicator function CT <= K
    VT = zeros(NPeriods*NReps, 1);  % CAT bond price (outputs)
    
    % placeholders for reusing densities in LR methods
    numer = zeros(NPeriods*NReps, 1);   % numerators
    denomILR = zeros(NPeriods*NReps, 1);    % denomimators for ILR
    denomMLR = zeros(NPeriods*NReps, 1);    % denominators for MLR
    
    % auxiliary variables for WLR (and ILR) implementation
    tempILR = zeros(NPeriods*NReps, 1);    % individual ILR values
    meanILR = zeros(NPeriods, 1);     % ILR estimate in each period
    weightWLR = zeros(NPeriods, 1);   % optimal weights for WLR
    
    % auxiliary variables for stochastic kriging
    Yhat = zeros(NPeriods, 1);   % estimated mean of simulation output
    Vhat = zeros(NPeriods, 1);   % estimated intrisic variance
    DsnX = zeros(NPeriods, 2); % place holder for design points
    
    for period = 1:NPeriods
        % index for current period
        ind_cur = (path-1)*NPeriods + period;
        % range of outputs for current experiment
        rng_cur = ((period-1)*NReps + 1):(period*NReps);
        % range of outputs for upto current experiment
        rng_cur_all = 1:period*NReps;
        % range of historical outputs
        rng_hist_all = 1:(period-1)*NReps;
        
        % current parameters
        lambda_cur = Paths(ind_cur, 1); theta_cur = Paths(ind_cur, 2);
        
        % computations of results (and storing results)
        NT(rng_cur) = poissrnd(lambda_cur, NReps, 1);
        CT(rng_cur) = gamrnd(NT(rng_cur), theta_cur, NReps, 1);
        IT(rng_cur) = (CT(rng_cur) <= K);
        VT(rng_cur) = p + (1-p) .* IT(rng_cur);
        
        % ===== SMC estimator
        PriceSMC(ind_cur) = mean(VT(rng_cur));
        
        % ===== Stochastic Kriging
        DsnX(period,:) = Paths(ind_cur, :);
        Yhat(period) = PriceSMC(ind_cur);
        Vhat(period) = var(VT(rng_cur))/NReps;
        B = [ones(period,1) ,DsnX(1:period,:)];
        
        if (period >= SKStart)
            if(period == SKStart) % cold start
                skriging_model = SK_GIBF_Fit(DsnX(1:period,:),Yhat(1:period),B,Vhat(1:period),[1 1]);
            else % warm start
                skriging_model = SK_GIBF_Fit(DsnX(1:period,:),Yhat(1:period),B,Vhat(1:period),[1,1],...
                    skriging_model.params);
            end
            
            % ===== GSK, estimate using "current" model
            PriceGSK(ind_cur) = SK_GIBF_Predict(skriging_model,DsnX(period,:),[1,DsnX(period,:)]);
            PriceGSK_mi(ind_cur) = SK_GIBF_Predict(skriging_model,[lambda_mi, theta_mi],[1, lambda_mi, theta_mi]);
            PriceGSK_hi(ind_cur) = SK_GIBF_Predict(skriging_model,[lambda_hi, theta_hi],[1, lambda_hi, theta_hi]);
        end
        
        % ===== Numerator desity, used in all LR methods
        numer(rng_cur_all) = besselpdf(CT(rng_cur_all), lambda_cur, theta_cur);
        
        % ===== ILR estimator
        denomILR(rng_cur) = numer(rng_cur);
        tempILR(rng_cur_all) = VT(rng_cur_all).*numer(rng_cur_all)./denomILR(rng_cur_all);
        for hist = 1:period
            rng_hist = ((hist-1)*NReps + 1):(hist*NReps);
            meanILR(hist) = mean(tempILR(rng_hist));
            weightWLR(hist) = 1/var(tempILR(rng_hist));
        end
        PriceILR(ind_cur) = mean(meanILR(1:period));
        
        % ==== WLR estimator (weighted ILR estimator)
        weightWLR(1:period) = weightWLR(1:period)/sum(weightWLR(1:period));   % optimal weights
        PriceWLR(ind_cur) = sum(meanILR(1:period).*weightWLR(1:period)); % WLR estimates
        
        % ===== MLR estimator
        for hist = 1:(period-1)
            ind_hist = hist + (path-1)*NPeriods;
            lambda_hist = Paths(ind_hist, 1);
            theta_hist = Paths(ind_hist, 2);
            
            denomMLR(rng_cur) = denomMLR(rng_cur) + besselpdf(CT(rng_cur), lambda_hist, theta_hist);
        end
        denomMLR(rng_cur_all) = denomMLR(rng_cur_all) + numer(rng_cur_all);
        PriceMLR(ind_cur) = period*mean(VT(rng_cur_all).*numer(rng_cur_all)./denomMLR(rng_cur_all));
        
        % ===== ILR  and WLR estimators for mi and hi states
        % denominator readily available, only need to recalculate numerator
        numer(rng_cur_all) = besselpdf(CT(rng_cur_all), lambda_mi, theta_mi);
        tempILR(rng_cur_all) = VT(rng_cur_all).*numer(rng_cur_all)./denomILR(rng_cur_all);
        for hist = 1:period
            rng_hist = ((hist-1)*NReps + 1):(hist*NReps);
            meanILR(hist) = mean(tempILR(rng_hist));
            weightWLR(hist) = 1/var(tempILR(rng_hist));
        end
        PriceILR_mi(ind_cur) = mean(meanILR(1:period));
        
        weightWLR(1:period) = weightWLR(1:period)/sum(weightWLR(1:period));   % optimal weights
        PriceWLR_mi(ind_cur) = sum(meanILR(1:period).*weightWLR(1:period)); % WLR estimates
        
        
        numer(rng_cur_all) = besselpdf(CT(rng_cur_all), lambda_hi, theta_hi);
        tempILR(rng_cur_all) = VT(rng_cur_all).*numer(rng_cur_all)./denomILR(rng_cur_all);
        for hist = 1:period
            rng_hist = ((hist-1)*NReps + 1):(hist*NReps);
            meanILR(hist) = mean(tempILR(rng_hist));
            weightWLR(hist) = 1/var(tempILR(rng_hist));
        end
        PriceILR_hi(ind_cur) = mean(meanILR(1:period));
        
        weightWLR(1:period) = weightWLR(1:period)/sum(weightWLR(1:period));   % optimal weights
        PriceWLR_hi(ind_cur) = sum(meanILR(1:period).*weightWLR(1:period)); % WLR estimates
        
        % ===== MLR estimators for mi and hi states
        numer(rng_cur_all) = besselpdf(CT(rng_cur_all), lambda_mi, theta_mi);
        PriceMLR_mi(ind_cur) = period*mean(VT(rng_cur_all).*numer(rng_cur_all)./denomMLR(rng_cur_all));
        
        numer(rng_cur_all) = besselpdf(CT(rng_cur_all), lambda_hi, theta_hi);
        PriceMLR_hi(ind_cur) = period*mean(VT(rng_cur_all).*numer(rng_cur_all)./denomMLR(rng_cur_all));
    end
    path
    toc
end

%% Estimate mean squared errors for different estimators
RMSE_SMC = reshape((PriceSMC - PriceTrue).^2, NPeriods, NPaths);
RMSE_ILR = reshape((PriceILR - PriceTrue).^2, NPeriods, NPaths);
RMSE_WLR = reshape((PriceWLR - PriceTrue).^2, NPeriods, NPaths);
RMSE_MLR = reshape((PriceMLR - PriceTrue).^2, NPeriods, NPaths);
RMSE_GSK = reshape((PriceGSK - PriceTrue).^2, NPeriods, NPaths);

RMSE_Y_SMC = mean(RMSE_SMC,2);  RMSE_E_SMC = std(RMSE_SMC, 0, 2)./sqrt(NPaths);
RMSE_Y_ILR = mean(RMSE_ILR,2);  RMSE_E_ILR = std(RMSE_ILR, 0, 2)./sqrt(NPaths);
RMSE_Y_WLR = mean(RMSE_WLR,2);  RMSE_E_WLR = std(RMSE_WLR, 0, 2)./sqrt(NPaths);
RMSE_Y_MLR = mean(RMSE_MLR,2);  RMSE_E_MLR = std(RMSE_MLR, 0, 2)./sqrt(NPaths);
RMSE_Y_GSK = mean(RMSE_GSK,2);  RMSE_E_GSK = std(RMSE_GSK, 0, 2)./sqrt(NPaths);

RMSE_ILR_mi = reshape((PriceILR_mi - Price_mi).^2, NPeriods, NPaths);
RMSE_ILR_hi = reshape((PriceILR_hi - Price_hi).^2, NPeriods, NPaths);
RMSE_Y_ILR_mi = mean(RMSE_ILR_mi,2);  RMSE_E_ILR_mi = std(RMSE_ILR_mi, 0, 2)./sqrt(NPaths);
RMSE_Y_ILR_hi = mean(RMSE_ILR_hi,2);  RMSE_E_ILR_hi = std(RMSE_ILR_hi, 0, 2)./sqrt(NPaths);

RMSE_WLR_mi = reshape((PriceWLR_mi - Price_mi).^2, NPeriods, NPaths);
RMSE_WLR_hi = reshape((PriceWLR_hi - Price_hi).^2, NPeriods, NPaths);
RMSE_Y_WLR_mi = mean(RMSE_WLR_mi,2);  RMSE_E_WLR_mi = std(RMSE_WLR_mi, 0, 2)./sqrt(NPaths);
RMSE_Y_WLR_hi = mean(RMSE_WLR_hi,2);  RMSE_E_WLR_hi = std(RMSE_WLR_hi, 0, 2)./sqrt(NPaths);

RMSE_MLR_mi = reshape((PriceMLR_mi - Price_mi).^2, NPeriods, NPaths);
RMSE_MLR_hi = reshape((PriceMLR_hi - Price_hi).^2, NPeriods, NPaths);
RMSE_Y_MLR_mi = mean(RMSE_MLR_mi,2);  RMSE_E_MLR_mi = std(RMSE_MLR_mi, 0, 2)./sqrt(NPaths);
RMSE_Y_MLR_hi = mean(RMSE_MLR_hi,2);  RMSE_E_MLR_hi = std(RMSE_MLR_hi, 0, 2)./sqrt(NPaths);

RMSE_GSK_mi = reshape((PriceGSK_mi - Price_mi).^2, NPeriods, NPaths);
RMSE_GSK_hi = reshape((PriceGSK_hi - Price_hi).^2, NPeriods, NPaths);
RMSE_Y_GSK_mi = mean(RMSE_GSK_mi,2);  RMSE_E_GSK_mi = std(RMSE_GSK_mi, 0, 2)./sqrt(NPaths);
RMSE_Y_GSK_hi = mean(RMSE_GSK_hi,2);  RMSE_E_GSK_hi = std(RMSE_GSK_hi, 0, 2)./sqrt(NPaths);

RMSE_Y_SMC_mi = CatBondPriceVar(lambda_mi, theta_mi, K, p, 40)/NReps*ones(NPeriods,1);
RMSE_Y_SMC_hi = CatBondPriceVar(lambda_hi, theta_hi, K, p, 40)/NReps*ones(NPeriods,1);

%% Plot results
% ----- Histogram of Parameters
figure
ha = tight_subplot(1,2,[.1 .05],[.1 .05],[.05 .05]);

axes(ha(1));
% histogram(Paths(:,1),200,'EdgeAlpha',0,'Normalization','pdf')
histogram(Paths(:,1),200,'Normalization','pdf')

axes(ha(2));
% histogram(Paths(:,2),200,'EdgeAlpha',0,'Normalization','pdf')
histogram(Paths(:,2),200,'Normalization','pdf')

% ----- Current state plot
figure
hold on

ref_x = logspace(0,log(100)/log(10),20);
ref_y =(ref_x.^(-1)*RMSE_Y_SMC(1));

plot(1:NPeriods, RMSE_Y_SMC, '-', 'color',[0,204,204]/255, 'LineWidth',4)
plot(1:NPeriods, RMSE_Y_ILR, '-.', 'color',[255,102,102]/255, 'LineWidth',4)
plot(1:NPeriods, RMSE_Y_WLR, '-.', 'color',[255,204,0]/255, 'LineWidth',4)
plot(1:NPeriods, RMSE_Y_MLR, '-.', 'color',[204,51,255]/255, 'LineWidth',4)
plot(SKStart:NPeriods,RMSE_Y_GSK(SKStart:NPeriods),':', 'color',[0,0,204]/255, 'LineWidth',4)
plot(ref_x,ref_y, '-', 'LineWidth',4, 'color',[0,0,0]/255)

% xlim([0,100])
% ylim([1e-2,2e1])

legend({'SMC','ILR','EWLR','MLR','GSK','Slope -1'}, 'FontSize',16, 'Location','southwest','FontWeight','bold')
legend boxoff
set(gca,'FontSize',14,'FontWeight','bold','yscale','log','xscale','log','FontName','FixedWidth','xscale','log','yscale','log')
xlabel('Number of periods, n', 'FontSize',15,'FontWeight','bold')
title('period n MSE for the current state','Position',[10 1e-3],'FontSize',15,'FontWeight','bold')
hold off

% ----- Fixed state plot
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
plot(SKStart:NPeriods,RMSE_Y_GSK_mi(SKStart:end), ':', 'color',[0,0,204]/255, 'LineWidth',4)
plot(ref_x,ref_y, '-', 'LineWidth',4, 'color',[0,0,0]/255)
% plot(2:NPeriods, RMSE_Y_Shr_mi(2:NPeriods), '-', 'color',[0.6,0.2,0], 'LineWidth',4)

set(gca,'FontSize',14,'FontWeight','bold','FontName','FixedWidth', ...
    'yscale','log','xscale','log', ...
    'YTickLabelMode','auto','XTickLabelMode','auto')
legend({'SMC','ILR','EWLR','MLR','GSK', 'Slope -1'}, 'FontSize',16, 'Location','southwest','FontWeight','bold')
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
plot(SKStart:NPeriods, RMSE_Y_GSK_hi(SKStart:end), ':', 'color',[0,0,204]/255, 'LineWidth',4)
plot(ref_x,ref_y, '-', 'LineWidth',4, 'color',[0,0,0]/255)

set(gca,'FontSize',14,'FontWeight','bold','FontName','FixedWidth', ...
    'yscale','log','xscale','log', ...
    'YTickLabelMode','auto','XTickLabelMode','auto')
legend({'SMC','ILR','EWLR','MLR','GSK', 'Slope -1'}, 'FontSize',16, 'Location','southwest','FontWeight','bold')
legend boxoff
xlabel('Number of periods, n', 'FontSize',15,'FontWeight','bold')
title('period n MSE for extreme state-x_{hi}','Position',[10 2e-2],'FontSize',15,'FontWeight','bold')
xlim([0,100])
ylim([3e-6,2e-2])
hold off