function Parameters = catbond_states(Model, NumPeriods, Stream)
% States of one sample path of the CAT bond example, as loss-model
% parameters (lambda_n, theta_n) = psi(X_n), n = 1, ..., N (Section 6.1)
% Input(s):
% -- Model: struct from catbond_settings.m
% -- NumPeriods: N, number of repeated experiments
% -- Stream: RandStream
% Output(s):
% -- Parameters: N x 2, row n = (lambda_n, theta_n)
% formula: X_1 ~ stationary distribution N(mu_inf/(1-phi), sigma^2/(1-phi^2));
%          X_n = mu_inf + phi X_{n-1} + eps_n;
%          psi(x) = lower + (upper - lower) / (1 + exp(-x)), per coordinate

StationaryMean = Model.LongRunMean ./ (1 - Model.AutoRegression);
StationaryStd = Model.NoiseStd ./ sqrt(1 - Model.AutoRegression.^2);
States = zeros(NumPeriods, 2);                                   % X_n
States(1, :) = StationaryMean + StationaryStd .* randn(Stream, 1, 2);
for n = 2:NumPeriods
    Noise = Model.NoiseStd .* randn(Stream, 1, 2);               % eps_n
    States(n, :) = Model.LongRunMean + Model.AutoRegression .* States(n - 1, :) + Noise;
end
Parameters = catbond_parameters(Model, States);
end
