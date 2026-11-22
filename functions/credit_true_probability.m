function Probability = credit_true_probability(Model, States)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% credit_true_probability.m                                                   %
% True value mu(x) = Pr(L_1 + L_2 > kappa | S_0 = x) of the credit risk       %
% example (Section 6.2), in closed form                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input(s):
% -- Model: struct from credit_settings.m
% -- States: N x 2 current asset values x
% Output(s):
% -- Probability: N x 1
% Derivation: the cost of a solvent company i is below its maximum
% a_i e^{-r(T-t)} Phi(-d_i(D_i)), so with these parameters a loss above kappa
% needs one company to default (cost a_j) and the other's cost to exceed
% kappa - a_j. Company i's cost exceeds c exactly when S_t,i < b_i(c), the
% asset value at which a_i e^{-r(T-t)} Phi(-d_i) = c (or b_i = D_i if c
% exceeds the maximum cost of a solvent company). The event is
%   {S_1 < D_1, S_2 < b_2} union {S_1 < b_1, S_2 < D_2},  b_i >= D_i,
% whose probability is, with C the t-copula cdf and F_i the lognormal cdfs,
%   C(F_1(D_1), F_2(b_2)) + C(F_1(b_1), F_2(D_2)) - C(F_1(D_1), F_2(D_2)).

Horizon = Model.Evaluation;                                   % T - t
Correlation = [1, Model.CopulaRho; Model.CopulaRho, 1];
a = Model.LossGivenDefault;
Bound = zeros(1, 2);                                          % b_1, b_2
for i = 1:2
    j = 3 - i;                                                % the defaulting company
    Cost = Model.Threshold - a(j);                            % c = kappa - a_j
    Quantile = Cost / (a(i) * exp(-Model.RiskFree * Horizon));   % Phi(-d_i) = Quantile
    if Quantile >= 1
        Bound(i) = 0;                                         % impossible
    else
        D2 = -norminv(Quantile);                              % d_i at the bound
        Bound(i) = Model.Debt(i) * exp(D2 * Model.Volatility(i) * sqrt(Horizon) ...
            - (Model.RiskFree - Model.Volatility(i)^2 / 2) * Horizon);
        Bound(i) = max(Bound(i), Model.Debt(i));              % default always counts
    end
end
% Check the premise: two solvent companies cannot exceed kappa together
% (the maximum cost of a solvent company is at S_t,i = D_i)
D2AtDebt = (Model.RiskFree - Model.Volatility.^2 / 2) * Horizon ./ (Model.Volatility * sqrt(Horizon));
MaxSolventCost = a .* exp(-Model.RiskFree * Horizon) .* normcdf(-D2AtDebt);
if sum(MaxSolventCost) > Model.Threshold
    error('credit_true_probability: closed form needs two solvent companies below kappa.');
end

NumStates = size(States, 1);
Probability = zeros(NumStates, 1);
for s = 1:NumStates
    [LogMean, LogStd] = scenario_marginals(Model, States(s, :));
    F = @(Value, i) logncdf(Value, LogMean(i), LogStd(i));    % marginal cdfs
    Both = copulacdf('t', [F(Model.Debt(1), 1), F(Model.Debt(2), 2)], Correlation, Model.CopulaDF);
    First = copulacdf('t', [F(Model.Debt(1), 1), F(Bound(2), 2)], Correlation, Model.CopulaDF);
    Second = copulacdf('t', [F(Bound(1), 1), F(Model.Debt(2), 2)], Correlation, Model.CopulaDF);
    Probability(s) = First + Second - Both;
end
end
