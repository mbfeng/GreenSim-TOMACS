function Prediction = sk_predict(Fit, Points)
% Stochastic kriging prediction at new points (Section 3, Equation 19)
% Input(s):
% -- Fit: struct from sk_fit.m
% -- Points: m x d prediction points (in the original units)
% Output(s):
% -- Prediction: m x 1, b(x)' beta + tau^2 r(x)' Sigma^{-1} (Y - B beta),
%    r(x) = correlations between x and the design points
Scaled = (Points - Fit.Lower) ./ Fit.Range;
NumPoints = size(Points, 1);
Prediction = zeros(NumPoints, 1);
for i = 1:NumPoints
    SquaredDistances = (Fit.Scaled - Scaled(i, :)).^2;              % k x d
    Correlation = exp(-SquaredDistances * Fit.Theta(:));            % r(x), k x 1
    Prediction(i) = [1, Points(i, :)] * Fit.Beta + Fit.Tau2 * (Correlation' * Fit.Weights);
end
end
