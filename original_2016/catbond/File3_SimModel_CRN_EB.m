% File3_SimModelCRN_Common.m
% Simulation model for CAT bond example
% Oct 1, 2016
% -- Implements four estimators using common random numbers (CRN)
%    1. Individual Likelihood Ratio (ILR)
%    2. Weighted ILR (WLR)
%    3. Multiple Likelihood Ratio (MLR)
%    4, Stochastic Kriging (STK)
clear all   % clear all variable in workspace
close all
rng(0);
clc
addpath([pwd, '\StochasticKriging'])
addpath([pwd, '\GIBF'])

% clearvars -except K p FaceVal ...
%     NPaths NPeriods NReps NTerms ...
%     Paths PriceTrue Price_lo Price_mi Price_hi ...
%     lambda_lo theta_lo lambda_hi theta_hi lambda_mi theta_mi;

File1_Initialization;
% File2_SimPaths;
%% Initialize placeholders for results
PriceSMC = zeros(NPeriods*NPaths, 1);   % Standard Monte Carlo estimates
PriceILR = zeros(NPeriods*NPaths, 1);   % ILR estimates
PriceWLR = zeros(NPeriods*NPaths, 1);   % Weighted ILR estimates
PriceMLR = zeros(NPeriods*NPaths, 1);   % MLR estimates
PriceSK1 = zeros(NPeriods*NPaths, 1);    % Stochastic Kriging estimates (current)
PriceSK2 = zeros(NPeriods*NPaths, 1);    % Stochastic Kriging estimates (predictive)
PriceShr = zeros(NPeriods*NPaths, 1);    % Shrinkage estimator

PriceILR_lo = zeros(NPeriods*NPaths, 1);    % ILR estimate for lo-params
PriceILR_mi = zeros(NPeriods*NPaths, 1);    % ILR estimate for mi-params
PriceILR_hi = zeros(NPeriods*NPaths, 1);    % ILR estimate for hi-params

PriceWLR_lo = zeros(NPeriods*NPaths, 1);    % WLR estimate for lo-params
PriceWLR_mi = zeros(NPeriods*NPaths, 1);    % WLR estimate for mi-params
PriceWLR_hi = zeros(NPeriods*NPaths, 1);    % WLR estimate for hi-params

PriceMLR_lo = zeros(NPeriods*NPaths, 1);    % MLR estimate for lo-params
PriceMLR_mi = zeros(NPeriods*NPaths, 1);    % MLR estimate for mi-params
PriceMLR_hi = zeros(NPeriods*NPaths, 1);    % MLR estimate for hi-params

PriceSK1_lo = zeros(NPeriods*NPaths, 1);    % SK1 estimate for lo-params
PriceSK1_mi = zeros(NPeriods*NPaths, 1);    % SK1 estimate for mi-params
PriceSK1_hi = zeros(NPeriods*NPaths, 1);    % SK1 estimate for hi-params

PriceShr_lo = zeros(NPeriods*NPaths, 1);    % Shrinkage estimate for lo-params
PriceShr_mi = zeros(NPeriods*NPaths, 1);    % Shrinkage estimate for mi-params
PriceShr_hi = zeros(NPeriods*NPaths, 1);    % Shrinkage estimate for hi-params

load('results_with_GIBF.mat')
rng(0);

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
    
    TempSMC = zeros(NPeriods, 1);   % estimated intrisic variance
    
    TempV_lo = zeros(NReps, 1);   % estimated intrisic variance
    TempV_mi = zeros(NReps, 1);   % estimated intrisic variance
    TempV_hi = zeros(NReps, 1);   % estimated intrisic variance
    
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
        
        % ===== Shrinkage estimator
        TempSMC(period) = PriceSMC(ind_cur);
        TempV_lo = p + (1-p) .* (gamrnd(poissrnd(lambda_lo, NReps, 1), theta_lo, NReps, 1) <= K);
        TempV_mi = p + (1-p) .* (gamrnd(poissrnd(lambda_mi, NReps, 1), theta_mi, NReps, 1) <= K);
        TempV_hi = p + (1-p) .* (gamrnd(poissrnd(lambda_hi, NReps, 1), theta_hi, NReps, 1) <= K);
        
        if (period == 1)
            PriceShr(ind_cur) = PriceSMC(ind_cur);
        else
            mu_est = mean(TempSMC(1:period));
            tau_est = var(TempSMC(1:period));
            
            sig_est = var(VT(rng_cur))/NReps;
            B_shrink = sig_est / (sig_est + tau_est);
            PriceShr(ind_cur) = B_shrink*mu_est + (1-B_shrink)*PriceSMC(ind_cur);
            
            sig_est = var(TempV_lo)/NReps;
            B_shrink = sig_est / (sig_est + tau_est);
            PriceShr_lo(ind_cur) = B_shrink*mu_est + (1-B_shrink)*mean(TempV_lo);
            
            sig_est = var(TempV_mi)/NReps;
            B_shrink = sig_est / (sig_est + tau_est);
            PriceShr_mi(ind_cur) = B_shrink*mu_est + (1-B_shrink)*mean(TempV_mi);
            
            sig_est = var(TempV_hi)/NReps;
            B_shrink = sig_est / (sig_est + tau_est);
            PriceShr_hi(ind_cur) = B_shrink*mu_est + (1-B_shrink)*mean(TempV_hi);
        end
        
%                 % ===== Stochastic Kriging
%                 DsnX(period,:) = Paths(ind_cur, :);
%                 Yhat(period) = PriceSMC(ind_cur);
%                 Vhat(period) = var(VT(rng_cur))/NReps;
%                 B = [ones(period,1) ,DsnX(1:period,:)];
%         
%                 if (period >= SKStart)
%                     if(period == SKStart) % cold start
%                         %                         skriging_model = SKfit(DsnX(1:period,:),Yhat(1:period),B,Vhat(1:period),2);
%                         skriging_model = SK_GIBF_Fit(DsnX(1:period,:),Yhat(1:period),B, Vhat(1:period), [1,1]);
%                     else % warm start
%                         %                         skriging_model = SKfit(DsnX(1:period,:),Yhat(1:period),B,Vhat(1:period),2,...
%                         %                             skriging_model.tausquared,skriging_model.theta);
%                         skriging_model = SK_GIBF_Fit(DsnX(1:period,:),Yhat(1:period),B, Vhat(1:period),[1,1],skriging_model.params);
%                     end
%         
%                     % ===== SK1, estimate using "current" model
%                     %                     PriceSK1(ind_cur) = SKpredict(skriging_model,DsnX(period,:),[1,DsnX(period,:)]);
%                     %                     PriceSK1_lo(ind_cur) = SKpredict(skriging_model,[lambda_lo, theta_lo],[1, lambda_lo, theta_lo]);
%                     %                     PriceSK1_mi(ind_cur) = SKpredict(skriging_model,[lambda_mi, theta_mi],[1, lambda_mi, theta_mi]);
%                     %                     PriceSK1_hi(ind_cur) = SKpredict(skriging_model,[lambda_hi, theta_hi],[1, lambda_hi, theta_hi]);
%                     PriceSK1(ind_cur) = SK_GIBF_Predict(skriging_model,DsnX(period,:),[1,DsnX(period,:)]);
%                     PriceSK1_lo(ind_cur) = SK_GIBF_Predict(skriging_model,[lambda_lo, theta_lo],[1, lambda_lo, theta_lo]);
%                     PriceSK1_mi(ind_cur) = SK_GIBF_Predict(skriging_model,[lambda_mi, theta_mi],[1, lambda_mi, theta_mi]);
%                     PriceSK1_hi(ind_cur) = SK_GIBF_Predict(skriging_model,[lambda_hi, theta_hi],[1, lambda_hi, theta_hi]);
%                     if(period < NPeriods)
%                         % ===== SK2, estimate using "previous" model
%                         %                         PriceSK2(ind_cur + 1) = SKpredict(skriging_model, Paths(ind_cur + 1, :), [1,Paths(ind_cur + 1, :)]);
%                         PriceSK2(ind_cur + 1) = SK_GIBF_Predict(skriging_model, Paths(ind_cur + 1, :), [1,Paths(ind_cur + 1, :)]);
%                     end
%                 end
        %
        %         % ===== Numerator desity, used in all LR methods
        %         %       Note: two different density functions
        %         %*** Compound pdf takes both number of claims and total claims
        %         %  numer(rng_cur_all) = compoundpdf(NT(rng_cur_all), CT(rng_cur_all), lambda_cur, theta_cur);
        %         %*** Bessel pdf takes only total claims
        %         numer(rng_cur_all) = besselpdf(CT(rng_cur_all), lambda_cur, theta_cur);
        %
        %         % ===== ILR estimator
        %         denomILR(rng_cur) = numer(rng_cur);
        %         tempILR(rng_cur_all) = VT(rng_cur_all).*numer(rng_cur_all)./denomILR(rng_cur_all);
        %         for hist = 1:period
        %             rng_hist = ((hist-1)*NReps + 1):(hist*NReps);
        %             meanILR(hist) = mean(tempILR(rng_hist));
        %             weightWLR(hist) = 1/var(tempILR(rng_hist));
        %         end
        %         PriceILR(ind_cur) = mean(meanILR(1:period));
        %
        %         % ==== WLR estimator (weighted ILR estimator)
        %         weightWLR(1:period) = weightWLR(1:period)/sum(weightWLR(1:period));   % optimal weights
        %         PriceWLR(ind_cur) = sum(meanILR(1:period).*weightWLR(1:period)); % WLR estimates
        %
        %         % ===== MLR estimator
        %         for hist = 1:(period-1)
        %             ind_hist = hist + (path-1)*NPeriods;
        %             lambda_hist = Paths(ind_hist, 1);
        %             theta_hist = Paths(ind_hist, 2);
        %
        %             %             denomMLR(rng_cur) = denomMLR(rng_cur) + compoundpdf(NT(rng_cur), CT(rng_cur), lambda_hist, theta_hist);
        %             denomMLR(rng_cur) = denomMLR(rng_cur) + besselpdf(CT(rng_cur), lambda_hist, theta_hist);
        %         end
        %         denomMLR(rng_cur_all) = denomMLR(rng_cur_all) + numer(rng_cur_all);
        %         PriceMLR(ind_cur) = period*mean(VT(rng_cur_all).*numer(rng_cur_all)./denomMLR(rng_cur_all));
        %
        %         % ===== ILR  and WLR estimators for lo, mi, and hi states
        %         % denominator readily available, only need to recalculate numerator
        %         %         numer(rng_cur_all) = compoundpdf(NT(rng_cur_all), CT(rng_cur_all), lambda_lo, theta_lo);
        %         numer(rng_cur_all) = besselpdf(CT(rng_cur_all), lambda_lo, theta_lo);
        %         tempILR(rng_cur_all) = VT(rng_cur_all).*numer(rng_cur_all)./denomILR(rng_cur_all);
        %         for hist = 1:period
        %             rng_hist = ((hist-1)*NReps + 1):(hist*NReps);
        %             meanILR(hist) = mean(tempILR(rng_hist));
        %             weightWLR(hist) = 1/var(tempILR(rng_hist));
        %         end
        %         PriceILR_lo(ind_cur) = mean(meanILR(1:period));
        %
        %         weightWLR(1:period) = weightWLR(1:period)/sum(weightWLR(1:period));   % optimal weights
        %         PriceWLR_lo(ind_cur) = sum(meanILR(1:period).*weightWLR(1:period)); % WLR estimates
        %
        %         %         numer(rng_cur_all) = compoundpdf(NT(rng_cur_all), CT(rng_cur_all), lambda_mi, theta_mi);
        %         numer(rng_cur_all) = besselpdf(CT(rng_cur_all), lambda_mi, theta_mi);
        %         tempILR(rng_cur_all) = VT(rng_cur_all).*numer(rng_cur_all)./denomILR(rng_cur_all);
        %         for hist = 1:period
        %             rng_hist = ((hist-1)*NReps + 1):(hist*NReps);
        %             meanILR(hist) = mean(tempILR(rng_hist));
        %             weightWLR(hist) = 1/var(tempILR(rng_hist));
        %         end
        %         PriceILR_mi(ind_cur) = mean(meanILR(1:period));
        %
        %         weightWLR(1:period) = weightWLR(1:period)/sum(weightWLR(1:period));   % optimal weights
        %         PriceWLR_mi(ind_cur) = sum(meanILR(1:period).*weightWLR(1:period)); % WLR estimates
        %
        %         %         numer(rng_cur_all) = compoundpdf(NT(rng_cur_all), CT(rng_cur_all), lambda_hi, theta_hi);
        %         numer(rng_cur_all) = besselpdf(CT(rng_cur_all), lambda_hi, theta_hi);
        %         tempILR(rng_cur_all) = VT(rng_cur_all).*numer(rng_cur_all)./denomILR(rng_cur_all);
        %         for hist = 1:period
        %             rng_hist = ((hist-1)*NReps + 1):(hist*NReps);
        %             meanILR(hist) = mean(tempILR(rng_hist));
        %             weightWLR(hist) = 1/var(tempILR(rng_hist));
        %         end
        %         PriceILR_hi(ind_cur) = mean(meanILR(1:period));
        %
        %         weightWLR(1:period) = weightWLR(1:period)/sum(weightWLR(1:period));   % optimal weights
        %         PriceWLR_hi(ind_cur) = sum(meanILR(1:period).*weightWLR(1:period)); % WLR estimates
        %
        %         % ===== MLR estimators for lo, mi, and hi states
        %         %         numer(rng_cur_all) = compoundpdf(NT(rng_cur_all), CT(rng_cur_all), lambda_lo, theta_lo);
        %         numer(rng_cur_all) = besselpdf(CT(rng_cur_all), lambda_lo, theta_lo);
        %         PriceMLR_lo(ind_cur) = period*mean(VT(rng_cur_all).*numer(rng_cur_all)./denomMLR(rng_cur_all));
        %
        %         %         numer(rng_cur_all) = compoundpdf(NT(rng_cur_all), CT(rng_cur_all), lambda_mi, theta_mi);
        %         numer(rng_cur_all) = besselpdf(CT(rng_cur_all), lambda_mi, theta_mi);
        %         PriceMLR_mi(ind_cur) = period*mean(VT(rng_cur_all).*numer(rng_cur_all)./denomMLR(rng_cur_all));
        %
        %         %         numer(rng_cur_all) = compoundpdf(NT(rng_cur_all), CT(rng_cur_all), lambda_hi, theta_hi);
        %         numer(rng_cur_all) = besselpdf(CT(rng_cur_all), lambda_hi, theta_hi);
        %         PriceMLR_hi(ind_cur) = period*mean(VT(rng_cur_all).*numer(rng_cur_all)./denomMLR(rng_cur_all));
    end
    path
    toc
end

%% Estimate mean squared errors for different estimators
RMSE_SMC = reshape((PriceSMC - PriceTrue(1:NPeriods*NPaths)).^2, NPeriods, NPaths);
RMSE_ILR = reshape((PriceILR - PriceTrue(1:NPeriods*NPaths)).^2, NPeriods, NPaths);
RMSE_WLR = reshape((PriceWLR - PriceTrue(1:NPeriods*NPaths)).^2, NPeriods, NPaths);
RMSE_MLR = reshape((PriceMLR - PriceTrue(1:NPeriods*NPaths)).^2, NPeriods, NPaths);
RMSE_SK1 = reshape((PriceSK1 - PriceTrue(1:NPeriods*NPaths)).^2, NPeriods, NPaths);
RMSE_SK2 = reshape((PriceSK2 - PriceTrue(1:NPeriods*NPaths)).^2, NPeriods, NPaths);
RMSE_Shr = reshape((PriceShr - PriceTrue(1:NPeriods*NPaths)).^2, NPeriods, NPaths);

RMSE_Y_SMC = mean(RMSE_SMC,2);  RMSE_E_SMC = std(RMSE_SMC, 0, 2)./sqrt(NPaths);
RMSE_Y_ILR = mean(RMSE_ILR,2);  RMSE_E_ILR = std(RMSE_ILR, 0, 2)./sqrt(NPaths);
RMSE_Y_WLR = mean(RMSE_WLR,2);  RMSE_E_WLR = std(RMSE_WLR, 0, 2)./sqrt(NPaths);
RMSE_Y_MLR = mean(RMSE_MLR,2);  RMSE_E_MLR = std(RMSE_MLR, 0, 2)./sqrt(NPaths);
RMSE_Y_SK1 = mean(RMSE_SK1,2);  RMSE_E_SK1 = std(RMSE_SK1, 0, 2)./sqrt(NPaths);
RMSE_Y_SK2 = mean(RMSE_SK2,2);  RMSE_E_SK2 = std(RMSE_SK2, 0, 2)./sqrt(NPaths);
RMSE_Y_Shr = mean(RMSE_Shr,2);  RMSE_E_Shr = std(RMSE_SK2, 0, 2)./sqrt(NPaths);

RMSE_ILR_lo = reshape((PriceILR_lo - Price_lo).^2, NPeriods, NPaths);
RMSE_ILR_mi = reshape((PriceILR_mi - Price_mi).^2, NPeriods, NPaths);
RMSE_ILR_hi = reshape((PriceILR_hi - Price_hi).^2, NPeriods, NPaths);
RMSE_Y_ILR_lo = mean(RMSE_ILR_lo,2);  RMSE_E_ILR_lo = std(RMSE_ILR_lo, 0, 2)./sqrt(NPaths);
RMSE_Y_ILR_mi = mean(RMSE_ILR_mi,2);  RMSE_E_ILR_mi = std(RMSE_ILR_mi, 0, 2)./sqrt(NPaths);
RMSE_Y_ILR_hi = mean(RMSE_ILR_hi,2);  RMSE_E_ILR_hi = std(RMSE_ILR_hi, 0, 2)./sqrt(NPaths);

RMSE_WLR_lo = reshape((PriceWLR_lo - Price_lo).^2, NPeriods, NPaths);
RMSE_WLR_mi = reshape((PriceWLR_mi - Price_mi).^2, NPeriods, NPaths);
RMSE_WLR_hi = reshape((PriceWLR_hi - Price_hi).^2, NPeriods, NPaths);
RMSE_Y_WLR_lo = mean(RMSE_WLR_lo,2);  RMSE_E_WLR_lo = std(RMSE_WLR_lo, 0, 2)./sqrt(NPaths);
RMSE_Y_WLR_mi = mean(RMSE_WLR_mi,2);  RMSE_E_WLR_mi = std(RMSE_WLR_mi, 0, 2)./sqrt(NPaths);
RMSE_Y_WLR_hi = mean(RMSE_WLR_hi,2);  RMSE_E_WLR_hi = std(RMSE_WLR_hi, 0, 2)./sqrt(NPaths);

RMSE_MLR_lo = reshape((PriceMLR_lo - Price_lo).^2, NPeriods, NPaths);
RMSE_MLR_mi = reshape((PriceMLR_mi - Price_mi).^2, NPeriods, NPaths);
RMSE_MLR_hi = reshape((PriceMLR_hi - Price_hi).^2, NPeriods, NPaths);
RMSE_Y_MLR_lo = mean(RMSE_MLR_lo,2);  RMSE_E_MLR_lo = std(RMSE_MLR_lo, 0, 2)./sqrt(NPaths);
RMSE_Y_MLR_mi = mean(RMSE_MLR_mi,2);  RMSE_E_MLR_mi = std(RMSE_MLR_mi, 0, 2)./sqrt(NPaths);
RMSE_Y_MLR_hi = mean(RMSE_MLR_hi,2);  RMSE_E_MLR_hi = std(RMSE_MLR_hi, 0, 2)./sqrt(NPaths);

RMSE_SK1_lo = reshape((PriceSK1_lo - Price_lo).^2, NPeriods, NPaths);
RMSE_SK1_mi = reshape((PriceSK1_mi - Price_mi).^2, NPeriods, NPaths);
RMSE_SK1_hi = reshape((PriceSK1_hi - Price_hi).^2, NPeriods, NPaths);
RMSE_Y_SK1_lo = mean(RMSE_SK1_lo,2);  RMSE_E_SK1_lo = std(RMSE_SK1_lo, 0, 2)./sqrt(NPaths);
RMSE_Y_SK1_mi = mean(RMSE_SK1_mi,2);  RMSE_E_SK1_mi = std(RMSE_SK1_mi, 0, 2)./sqrt(NPaths);
RMSE_Y_SK1_hi = mean(RMSE_SK1_hi,2);  RMSE_E_SK1_hi = std(RMSE_SK1_hi, 0, 2)./sqrt(NPaths);

RMSE_Shr_lo = reshape((PriceShr_lo - Price_lo).^2, NPeriods, NPaths);
RMSE_Shr_mi = reshape((PriceShr_mi - Price_mi).^2, NPeriods, NPaths);
RMSE_Shr_hi = reshape((PriceShr_hi - Price_hi).^2, NPeriods, NPaths);
RMSE_Y_Shr_lo = mean(RMSE_Shr_lo,2);  RMSE_E_Shr_lo = std(RMSE_Shr_lo, 0, 2)./sqrt(NPaths);
RMSE_Y_Shr_mi = mean(RMSE_Shr_mi,2);  RMSE_E_Shr_mi = std(RMSE_Shr_mi, 0, 2)./sqrt(NPaths);
RMSE_Y_Shr_hi = mean(RMSE_Shr_hi,2);  RMSE_E_Shr_hi = std(RMSE_Shr_hi, 0, 2)./sqrt(NPaths);

RMSE_Y_SMC_lo = CatBondPriceVar(lambda_lo, theta_lo, K, p, 40)/NReps*ones(NPeriods,1);
RMSE_Y_SMC_mi = CatBondPriceVar(lambda_mi, theta_mi, K, p, 40)/NReps*ones(NPeriods,1);
RMSE_Y_SMC_hi = CatBondPriceVar(lambda_hi, theta_hi, K, p, 40)/NReps*ones(NPeriods,1);

save('tempResults.mat')