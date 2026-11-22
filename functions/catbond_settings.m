function Model = catbond_settings()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% catbond_settings.m                                                          %
% Parameters of the catastrophe bond example (Section 6.1 of Feng and Staum, %
% ACM TOMACS 2017)                                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output(s): struct Model with
% -- .Trigger: K, trigger level of the insured loss (thousands of millions)
% -- .Recovery: p, fraction of face value paid if the loss exceeds K
% -- .LongRunMean, .AutoRegression, .NoiseStd: AR(1) state process
%    X_n = mu_inf + phi X_{n-1} + eps_n, eps_n ~ N(0, diag(NoiseStd.^2))
% -- .LambdaRange, .ThetaRange: bounds of the parameters (lambda, theta)
%    = psi(x) (expected number of hurricanes, expected loss per hurricane)
% -- .SeriesTerms: terms of the Poisson-Gamma series for the true price
% -- .FixedStates: struct with the parameters of the fixed states x_lo, x_mi,
%    x_hi (corner, centre and corner of the parameter range)

Model.Trigger = 25;                 % K
Model.Recovery = 0.5;               % p
Model.LongRunMean = [0, 0];         % mu_inf
Model.AutoRegression = [0.6, 0.5];  % phi
Model.NoiseStd = [0.8, 0.5];        % standard deviations of eps_n
Model.LambdaRange = [2, 4];         % (lambda_lower, lambda_upper)
Model.ThetaRange = [4, 6];          % (theta_lower, theta_upper)
Model.SeriesTerms = 40;             % m = 0, ..., 40 in the true price
Model.FixedStates.lo = [2, 4];      % parameters (lambda, theta) of x_lo
Model.FixedStates.mi = [3, 5];      % x_mi = (0, 0)
Model.FixedStates.hi = [4, 6];      % x_hi = (inf, inf), a limit point
end
