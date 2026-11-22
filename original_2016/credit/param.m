% param.m
% generate parameters in the current workspace

%% initializing workspace
clear;      % clear variables
close all;  % close all figures
clc;        % clear screen
rng(0);     % random number seed

%% ===== 1 KMV model parameters
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

%% Function that calculates true probability
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

%% ===== 2. Experiment design parameters
NAsset = 2;      % number of assets
NReps = 1e3;   % number of prediction pts at time t
NPaths = 1e4;    % number of sample paths
NPeriods = 52;  % number of periods in each sample path
SKStart = 5;    % starting period of stochastic kriging

clearvars -except Port TrueProbFunc ...
    NAsset NReps NPaths NPeriods SKStart
