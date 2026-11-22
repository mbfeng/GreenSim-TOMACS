function [Price, SMCVariance] = catbond_price(Model, Parameters, NumReplications)
% True price mu(x) of the CAT bond per dollar (Equation 26) and the variance
% of the standard Monte Carlo estimator with r replications
% Input(s):
% -- Model: struct from catbond_settings.m
% -- Parameters: N x 2, rows (lambda, theta)
% -- NumReplications: r (for the SMC variance; optional)
% Output(s):
% -- Price: N x 1, mu(x) = p + (1 - p) sum_{m>=0} Pois(m; lambda) G(K; m, theta)
%    (series truncated after Model.SeriesTerms terms; G = Gamma cdf, and
%    G(K; 0, theta) = 1)
% -- SMCVariance: N x 1, (1 - p)^2 P (1 - P) / r, P = Pr(Y <= K)

Terms = 0:Model.SeriesTerms;
NumStates = size(Parameters, 1);
Probability = zeros(NumStates, 1);                     % Pr(Y <= K)
for s = 1:NumStates
    Weights = poisspdf(Terms, Parameters(s, 1));                      % Pr(M = m)
    Conditional = gamcdf(Model.Trigger, Terms, Parameters(s, 2));     % Pr(Y <= K | M = m)
    Probability(s) = Weights * Conditional';
end
Price = Model.Recovery + (1 - Model.Recovery) * Probability;
if nargin > 2
    SMCVariance = (1 - Model.Recovery)^2 * Probability .* (1 - Probability) / NumReplications;
end
end
