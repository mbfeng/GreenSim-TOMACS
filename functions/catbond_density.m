function Density = catbond_density(Losses, Parameters)
% Density h(y; x) of the insured loss Y given the parameters (lambda, theta)
% of a state (Equation 27, Appendix A): a point mass at 0 and a
% continuous part on y > 0
% Input(s):
% -- Losses: column vector of losses y
% -- Parameters: 1 x 2, (lambda, theta)
% Output(s):
% -- Density: h(y; x); for y = 0 the probability e^{-lambda} of no claims
% formula (y > 0): h(y) = lambda I_1(2 sqrt(lambda y / theta)) /
%                  (theta e^{lambda + y/theta} sqrt(lambda y / theta)),
%   I_1 = modified Bessel function of the first kind of order 1
% The likelihood ratios h(y; x) / h(y; x') of all estimators use this
% density with respect to (point mass at 0) + (Lebesgue measure on y > 0).

Lambda = Parameters(1); Theta = Parameters(2);
Density = zeros(size(Losses));
Positive = Losses > 0;
Argument = sqrt(Lambda * Losses(Positive) / Theta);
Density(Positive) = Lambda * besseli(1, 2 * Argument) ...
    ./ (Theta * exp(Lambda + Losses(Positive) / Theta) .* Argument);
Density(~Positive) = exp(-Lambda);                 % Pr(M = 0)
end
