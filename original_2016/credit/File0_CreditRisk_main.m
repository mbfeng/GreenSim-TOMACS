% File0_CreditRisk_main.m
% Green simulation for credit risk example
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

%% Experiment Parameters
% ----- Experiment design parameters
NReps = 1e3;   % number of prediction pts at time t
NPaths = 1e4;    % number of sample paths
NPeriods = 52;  % number of periods in each sample path
SKStart = 5;    % starting period of stochastic kriging

% ----- Simulation model parameters
Port.T = 5;     % asset valuation period = 5 years
Port.t = 0.5;   % projection period = 6 months
Port.dt = 1/52; % re-evaluation period = 1 week

V0_Eqty = [15, 5];  % value of equity
V0_Debt = [85, 85]; % value of debt
V0_Firm = V0_Debt + V0_Eqty;   % firm values at time 0
Port.S0 = V0_Firm;
Port.K = V0_Debt;
Port.LGD = [5, 4];  % loss given default

Port.mu = [0.15, 0.10];     % drift
Port.rf = 0.05;    % risk free rate
Port.vol = [0.3, 0.2];     % volatility
Port.Threshold = 6;    % threshold for large loss

Port.copula = 't';      % choice of copula
Port.df = 3;            % degree of freedom in t-copula simulation
Port.Cor = [1     0.5;
    0.5   1;];  % correlation matrix

Port.mu_t = (Port.mu - (Port.vol.^ 2) ./ 2) .* Port.t;
Port.sig_t = Port.vol .* sqrt(Port.t);
Port.mu_delt = (Port.mu - (Port.vol.^ 2) ./ 2) .* Port.dt;
Port.sig_delt = Port.vol .* sqrt(Port.dt);

% ----- function that calculates true probability
options = optimoptions('fsolve','TolFun',1e-16);
d2 = @(S, K, vol) (log(S/K) + (Port.rf - Port.t*vol^2) *(Port.T-Port.t)) / (vol*sqrt(Port.T-Port.t));
LossFunc = @(x, K, vol, LGD) exp(-Port.rf * (Port.T-Port.t)) * normcdf(-d2(x,K,vol)) * LGD;

[Stock1Bound, ~, flag1, ~] = fsolve(@(x) LossFunc(x,Port.K(1), Port.vol(1), Port.LGD(1)) - (Port.Threshold - Port.LGD(2)), Port.S0(1), options);
if(flag1 < 0 || Stock1Bound<Port.K(1))   
    flag1 = 0;
else
    flag1 = 1;
end
[Stock2Bound, ~, flag2, ~] = fsolve(@(x) LossFunc(x,Port.K(2), Port.vol(2), Port.LGD(2)) - (Port.Threshold - Port.LGD(1)), Port.S0(2), options);
if(flag2 < 0 || Stock2Bound<Port.K(2))    
    flag2 = 0;
else
    flag2 = 1;
end

TrueProbFunc = @(Price1, Price2) flag2*copulacdf('t', [logncdf(Port.K(1), log(Price1) + Port.mu_t(1), Port.sig_t(1)), logncdf(Stock2Bound, log(Price2) + Port.mu_t(2), Port.sig_t(2))], Port.Cor,Port.df)...
    + flag1*copulacdf('t', [logncdf(Stock1Bound, log(Price1) + Port.mu_t(1), Port.sig_t(1)), logncdf(Port.K(2), log(Price2) + Port.mu_t(2), Port.sig_t(2))], Port.Cor,Port.df)...
    - flag1*flag2*copulacdf('t', [logncdf(Port.K(1), log(Price1) + Port.mu_t(1), Port.sig_t(1)), logncdf(Port.K(2), log(Price2) + Port.mu_t(2), Port.sig_t(2))], Port.Cor,Port.df);

%% Simulate sample paths
PathsMu = norminv(copularnd('t', Port.Cor, Port.df, NPeriods * NPaths));
TrueProb = zeros(NPeriods * NPaths, 1);

path = 0;
for ind = 1:(NPeriods * NPaths)
    if(mod(ind,NPeriods) == 1)
        PathsMu(ind,:) = log(Port.S0);
        path = path + 1
    else
        PathsMu(ind,:) = PathsMu(ind-1,:) + Port.mu_delt + Port.sig_delt.*PathsMu(ind,:);
    end
    TrueProb(ind) = TrueProbFunc(exp(PathsMu(ind,1)), exp(PathsMu(ind,2)));
end

%% Simulation model for credit risk evaluation
% ----- placeholders for results
PriceSMC = zeros(NPeriods*NPaths, 1);   % Standard Monte Carlo (SMC)
PriceILR = zeros(NPeriods*NPaths, 1);   % Individual Likelihood Ratio (ILR)
PriceWLR = zeros(NPeriods*NPaths, 1);   % Weighted Likelihood Ratio (WLR)
PriceMLR = zeros(NPeriods*NPaths, 1);   % Mixture Likelihood Ratio (MLR)

for path = 1:NPaths
    tic
    
    % The r.v's and the function value pairs
    Vt_Firm = zeros(NReps * NPeriods, 2);
    BigLossInd = zeros(NReps * NPeriods, 1);
    
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
        mu1_cur = PathsMu(ind_cur,1); mu2_cur = PathsMu(ind_cur,2);
        
        % Generate correlated uniform(0,1) using copula
        U_0t = copularnd('t', Port.Cor, Port.df, NReps);
        Vt_Firm(rng_cur,1) = logninv(U_0t(:,1), mu1_cur + Port.mu_t(1), Port.sig_t(1));
        Vt_Firm(rng_cur,2) = logninv(U_0t(:,2), mu2_cur + Port.mu_t(2), Port.sig_t(2));
        BigLossInd(rng_cur) = LossAtt(Vt_Firm(rng_cur,:), Port) > Port.Threshold;
        
        % ===== SMC estimator
        PriceSMC(ind_cur) = mean(BigLossInd(rng_cur));
        
        % ===== Numerator desity, used in all LR methods
        numer(rng_cur_all) = CopulaPDF(Vt_Firm(rng_cur_all,1), Vt_Firm(rng_cur_all,2), mu1_cur, mu2_cur, Port);
        numer(isnan(numer)) = 0;
        
        % ===== ILR estimator
        denomILR(rng_cur) = numer(rng_cur);
        tempILR(rng_cur_all) = BigLossInd(rng_cur_all).*numer(rng_cur_all)./denomILR(rng_cur_all);
        for hist = 1:period
            rng_hist = ((hist-1)*NReps + 1):(hist*NReps);
            meanILR(hist) = mean(tempILR(rng_hist));
            if (var(tempILR(rng_hist)) > 0)
                weightWLR(hist) = 1/var(tempILR(rng_hist));
            else
                weightWLR(hist) = 0;
            end
        end
        PriceILR(ind_cur) = mean(meanILR(1:period));
        
        % ==== WLR estimator (weighted ILR estimator)
        if (sum(weightWLR(1:period)) == 0)
            weightWLR(1:period) = 0;
        else
            weightWLR(1:period) = weightWLR(1:period)/sum(weightWLR(1:period));   % optimal weights
        end
        PriceWLR(ind_cur) = sum(meanILR(1:period).*weightWLR(1:period)); % WLR estimates
        
        % ===== MLR estimator
        for hist = 1:(period-1)
            ind_hist = hist + (path-1)*NPeriods;
            mu1_hist = PathsMu(ind_hist, 1);
            mu2_hist = PathsMu(ind_hist, 2);
            
            
            denomMLR(rng_cur) = denomMLR(rng_cur) + CopulaPDF(Vt_Firm(rng_cur,1), Vt_Firm(rng_cur,2), mu1_hist, mu2_hist,Port);
        end
        denomMLR(rng_cur_all) = denomMLR(rng_cur_all) + numer(rng_cur_all);
        PriceMLR(ind_cur) = period*mean(BigLossInd(rng_cur_all).*numer(rng_cur_all)./denomMLR(rng_cur_all));
    end
    path
    toc
end

RMSE_SMC = reshape((PriceSMC - TrueProb(1:NPeriods*NPaths)).^2, NPeriods, NPaths);
RMSE_ILR = reshape((PriceILR - TrueProb(1:NPeriods*NPaths)).^2, NPeriods, NPaths);
RMSE_WLR = reshape((PriceWLR - TrueProb(1:NPeriods*NPaths)).^2, NPeriods, NPaths);
RMSE_MLR = reshape((PriceMLR - TrueProb(1:NPeriods*NPaths)).^2, NPeriods, NPaths);

RMSE_Y_SMC = mean(RMSE_SMC,2);  RMSE_E_SMC = std(RMSE_SMC, 0, 2)./sqrt(NPaths);
RMSE_Y_ILR = mean(RMSE_ILR,2);  RMSE_E_ILR = std(RMSE_ILR, 0, 2)./sqrt(NPaths);
RMSE_Y_WLR = mean(RMSE_WLR,2);  RMSE_E_WLR = std(RMSE_WLR, 0, 2)./sqrt(NPaths);
RMSE_Y_MLR = mean(RMSE_MLR,2);  RMSE_E_MLR = std(RMSE_MLR, 0, 2)./sqrt(NPaths);

%% Plot results

figure
hold on

errorbar(RMSE_Y_SMC, RMSE_E_SMC, '-', 'color',[0,204,204]/255, 'LineWidth',3)
errorbar(RMSE_Y_ILR, RMSE_E_ILR, '-.', 'color',[255,102,102]/255, 'LineWidth',3)
errorbar(RMSE_Y_WLR, RMSE_E_WLR, '-.', 'color',[255,204,0]/255, 'LineWidth',3)
errorbar(RMSE_Y_MLR, RMSE_E_MLR, '-.', 'color',[204,51,255]/255, 'LineWidth',3)

set(gca,'FontSize',14,'FontWeight','bold','yscale','log','xscale','log','FontName','FixedWidth','xscale','log','yscale','log')
xlim([1,52])
ylim([1e-6,1e-3])
legend({'SMC','ILR','EWLR','MLR'}, 'FontSize',16, 'Location','southwest','FontWeight','bold')
legend boxoff
xlabel('Number of periods, n', 'FontSize',15,'FontWeight','bold')
title('period n variance for the current state','FontSize',15,'FontWeight','bold')
hold off

