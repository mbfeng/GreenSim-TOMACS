function Fit = sk_fit(Design, Responses, IntrinsicVariances, WarmStart)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sk_fit.m                                                                    %
% Stochastic kriging (Ankenman, Nelson and Staum 2010) with a linear trend    %
% and Gaussian correlation, fitted by maximum profile likelihood (Section 3)  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input(s):
% -- Design: k x d design points (the parameters of experiments 1, ..., n)
% -- Responses: k x 1 simulation outputs at the design points (SMC estimates)
% -- IntrinsicVariances: k x 1 variances of the responses (sample variance /r)
% -- WarmStart: optional struct with fields Tau2 and Theta of an earlier fit,
%    used as the starting point of the optimization
% Output(s): struct Fit with the estimates Tau2 (extrinsic variance tau^2),
% Theta (1 x d correlation parameters), Beta (trend coefficients) and what
% sk_predict.m needs (scaled design, Cholesky factor, residual solve)
% Model: Y(x) = b(x)' beta + M(x) + eps, b(x) = [1, x], M a Gaussian process
%   with covariance tau^2 exp(-sum_p theta_p (u_p - u'_p)^2) in the design
%   scaled to [0, 1]^d, eps independent with known variances V.
% formula: Sigma = tau^2 R(theta) + diag(V); beta = GLS estimate given Sigma;
%   (tau^2, theta) minimize the negative log profile likelihood
%   0.5 log det(Sigma) + 0.5 (Y - B beta)' Sigma^{-1} (Y - B beta),
%   subject to tau^2 >= 0.001 tau^2_0 and theta >= 0.001, where tau^2_0 is
%   the starting value: the variance of the ordinary least-squares
%   residuals (cold start) or the warm-start value.

[NumPoints, Dim] = size(Design);
Lower = min(Design, [], 1);
Range = max(Design, [], 1) - Lower;
Scaled = (Design - Lower) ./ Range;                     % design in [0, 1]^d
Basis = [ones(NumPoints, 1), Design];                   % b(x) = [1, x]

% Squared distances per dimension: SquaredDistances(:, :, p) = (u_ip - u_jp)^2
SquaredDistances = zeros(NumPoints, NumPoints, Dim);
for p = 1:Dim
    SquaredDistances(:, :, p) = (Scaled(:, p) - Scaled(:, p)').^2;
end

% ===== Starting values
if nargin < 4 || isempty(WarmStart)
    OLS = Basis \ Responses;                            % ordinary least squares
    Tau2Start = var(Responses - Basis * OLS);
    % correlation 1/2 at the average distance between distinct design points
    Upper = triu(true(NumPoints), 1);
    AverageDistance = zeros(1, Dim);
    for p = 1:Dim
        Distances = sqrt(SquaredDistances(:, :, p));
        AverageDistance(p) = mean(Distances(Upper));
    end
    ThetaStart = (log(2) / Dim) * AverageDistance.^(-2);
else
    Tau2Start = WarmStart.Tau2;
    ThetaStart = WarmStart.Theta;
end

% ===== Maximum profile likelihood
Objective = @(Parameters) negative_log_likelihood(Parameters, SquaredDistances, ...
    Basis, IntrinsicVariances, Responses);
LowerBounds = [0.001 * Tau2Start, 0.001 * ones(1, Dim)];
Options = optimoptions('fmincon', 'Display', 'off', 'MaxFunctionEvaluations', 1e6, ...
    'MaxIterations', 1000);
% fmincon's interior-point steps can be ill-conditioned near the bounds; its
% internal "nearly singular" warnings do not affect the optimum, so they are
% silenced here (and the previous setting restored)
PreviousWarning = warning('off', 'MATLAB:nearlySingularMatrix');
Best = fmincon(Objective, [Tau2Start, ThetaStart], [], [], [], [], LowerBounds, [], [], Options);
warning(PreviousWarning);

% ===== Estimates and quantities for prediction
Fit.Tau2 = Best(1);
Fit.Theta = Best(2:end);
[~, Fit.Beta, Fit.Cholesky, Fit.Weights] = negative_log_likelihood(Best, SquaredDistances, ...
    Basis, IntrinsicVariances, Responses);
Fit.Scaled = Scaled;
Fit.Lower = Lower;
Fit.Range = Range;
end


function [Value, Beta, Cholesky, Weights] = negative_log_likelihood(Parameters, SquaredDistances, Basis, IntrinsicVariances, Responses)
% Negative log profile likelihood of (tau^2, theta); also returns the GLS
% trend, the Cholesky factor of Sigma and Sigma^{-1} (Y - B beta)
Tau2 = Parameters(1);
Theta = reshape(Parameters(2:end), 1, 1, []);
Correlation = exp(-sum(SquaredDistances .* Theta, 3));          % R(theta)
Sigma = Tau2 * Correlation + diag(IntrinsicVariances);
[Cholesky, NotPositive] = chol(Sigma, 'lower');
if NotPositive
    Value = Inf; Beta = []; Weights = [];
    return
end
WhitenedBasis = Cholesky \ Basis;                               % L^{-1} B
WhitenedResponses = Cholesky \ Responses;                       % L^{-1} Y
Beta = WhitenedBasis \ WhitenedResponses;                        % GLS estimate
Whitened = WhitenedResponses - WhitenedBasis * Beta;            % L^{-1}(Y - B beta)
Value = sum(log(diag(Cholesky))) + 0.5 * (Whitened' * Whitened);
Weights = Cholesky' \ Whitened;                                 % Sigma^{-1}(Y - B beta)
end
