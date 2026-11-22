%% MS_IS2.m
% Implement the second type of importance sampling

addpath([pwd, '/StochasticKriging'])
warning('off' , 'all')
%% Load parameters
% param;
% genpaths;
rng(0);

%% Initialize placeholders for results
% PriceSMC = zeros(NPeriods*NPaths, 1);   % Standard Monte Carlo estimates
% PriceILR = zeros(NPeriods*NPaths, 1);   % ILR estimates
% PriceWLR = zeros(NPeriods*NPaths, 1);   % Weighted ILR estimates
% PriceMLR = zeros(NPeriods*NPaths, 1);   % MLR estimates
% PriceSK1 = zeros(NPeriods*NPaths, 1);    % Stochastic Kriging estimates (current)
% PriceSK2 = zeros(NPeriods*NPaths, 1);    % Stochastic Kriging estimates (predictive)
load('results_nova.mat')
PriceShr = zeros(NPeriods*NPaths, 1);    % Stochastic Kriging estimates (predictive)
for path = 1:NPaths
    tic
    
    % The r.v's and the function value pairs
    Vt_Firm = zeros(NReps * NPeriods, NAsset);
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
        
%         % ===== Stochastic Kriging
%         DsnX(period,:) = PathsMu(ind_cur, :);
%         Yhat(period) = 1e2*PriceSMC(ind_cur);
%         Vhat(period) = 1e4*var(BigLossInd(rng_cur))/NReps;
%         B = [ones(period,1) ,DsnX(1:period,:)];
%         
%         if (period >= SKStart)
%             if(period == SKStart) % cold start
%                 skriging_model = SKfit(DsnX(1:period,:),Yhat(1:period),B,Vhat(1:period),2);
%             else % warm start
%                 skriging_model = SKfit(DsnX(1:period,:),Yhat(1:period),B,Vhat(1:period),2,...
%                     skriging_model.tausquared,skriging_model.theta);
%             end
%             
%             % ===== SK1, estimate using "current" model
%             PriceSK1(ind_cur) = SKpredict(skriging_model,DsnX(period,:),[1,DsnX(period,:)]);
%             
%             if(period < NPeriods)
%                 % ===== SK2, estimate using "previous" model
%                 PriceSK2(ind_cur + 1) = SKpredict(skriging_model, PathsMu(ind_cur + 1, :), [1,PathsMu(ind_cur + 1, :)]);
%             end
%         end
        
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

%% Estimate mean squared errors for different estimators
RMSE_SMC = reshape((PriceSMC - TrueProb(1:NPeriods*NPaths)).^2, NPeriods, NPaths);
RMSE_ILR = reshape((PriceILR - TrueProb(1:NPeriods*NPaths)).^2, NPeriods, NPaths);
RMSE_WLR = reshape((PriceWLR - TrueProb(1:NPeriods*NPaths)).^2, NPeriods, NPaths);
RMSE_MLR = reshape((PriceMLR - TrueProb(1:NPeriods*NPaths)).^2, NPeriods, NPaths);
% RMSE_SK1 = reshape((PriceSK1 - PriceTrue(1:NPeriods*NPaths)).^2, NPeriods, NPaths);
% RMSE_SK2 = reshape((PriceSK2 - PriceTrue(1:NPeriods*NPaths)).^2, NPeriods, NPaths);

RMSE_Y_SMC = mean(RMSE_SMC,2);  RMSE_E_SMC = std(RMSE_SMC, 0, 2)./sqrt(NPaths);
RMSE_Y_ILR = mean(RMSE_ILR,2);  RMSE_E_ILR = std(RMSE_ILR, 0, 2)./sqrt(NPaths);
RMSE_Y_WLR = mean(RMSE_WLR,2);  RMSE_E_WLR = std(RMSE_WLR, 0, 2)./sqrt(NPaths);
RMSE_Y_MLR = mean(RMSE_MLR,2);  RMSE_E_MLR = std(RMSE_MLR, 0, 2)./sqrt(NPaths);
% RMSE_Y_SK1 = mean(RMSE_SK1,2);  RMSE_E_SK1 = std(RMSE_SK1, 0, 2)./sqrt(NPaths);
% RMSE_Y_SK2 = mean(RMSE_SK2,2);  RMSE_E_SK2 = std(RMSE_SK2, 0, 2)./sqrt(NPaths);


loglog(1:NPeriods, RMSE_Y_SMC/1e4,...
    1:NPeriods, RMSE_Y_ILR,...
    1:NPeriods, RMSE_Y_WLR,...
    1:NPeriods, RMSE_Y_MLR)

figure
hold on
errorbar(RMSE_Y_SMC,RMSE_E_SMC)
errorbar(RMSE_Y_ILR,RMSE_E_ILR)
errorbar(RMSE_Y_WLR,RMSE_E_WLR)
errorbar(RMSE_Y_MLR,RMSE_E_MLR)
set(gca,'yscale','log','xscale','log')

save('results_nova.mat')