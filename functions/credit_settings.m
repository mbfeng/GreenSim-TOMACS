function Model = credit_settings()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% credit_settings.m                                                           %
% Parameters of the credit risk example (Section 6.2 of Feng and Staum,       %
% ACM TOMACS 2017)                                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output(s): struct Model with
% -- .InitialAssets: 1 x 2, asset values (100, 90) at the first experiment
% -- .Debt: D; .LossGivenDefault: a = (5, 4)
% -- .Drift, .Volatility: eta = (15%, 10%), varsigma = (30%, 20%)
% -- .RiskFree: r = 5%
% -- .CopulaRho, .CopulaDF: t-copula with 3 degrees of freedom, correlation 0.5
% -- .RepeatPeriod: Delta t = 1/52 (weekly experiments)
% -- .Projection: t = 0.5, risk horizon
% -- .Evaluation: T - t = 4.5, time from the horizon to the loans' maturity
%    T = 5 (the digital put in the cost of an anticipated default)
% -- .Threshold: kappa = 6

Model.InitialAssets = [100, 90];
Model.Debt = [85, 85];
Model.LossGivenDefault = [5, 4];
Model.Drift = [0.15, 0.10];
Model.Volatility = [0.30, 0.20];
Model.RiskFree = 0.05;
Model.CopulaRho = 0.5;
Model.CopulaDF = 3;
Model.RepeatPeriod = 1/52;
Model.Projection = 0.5;
Model.Evaluation = 5 - 0.5;          % T - t
Model.Threshold = 6;
end
