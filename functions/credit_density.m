function Density = credit_density(Model, Scenarios, State)
% Density h(y; x) of the scenarios given a state: t-copula density times the
% lognormal marginal densities (used in the likelihood ratios, Section 2)
% Input(s):
% -- Model: struct from credit_settings.m
% -- Scenarios: N x 2 asset values y
% -- State: 1 x 2, state x
% Output(s):
% -- Density: N x 1, h(y; x)
% formula: h(y; x) = c(F_1(y_1), F_2(y_2)) f_1(y_1) f_2(y_2)

Correlation = [1, Model.CopulaRho; Model.CopulaRho, 1];
[LogMean, LogStd] = scenario_marginals(Model, State);
Uniforms = zeros(size(Scenarios));      % marginal CDF values F_i(y_i)
Marginals = ones(size(Scenarios, 1), 1); % product of marginal densities
for Asset = 1:2
    Uniforms(:, Asset) = logncdf(Scenarios(:, Asset), LogMean(Asset), LogStd(Asset));
    Marginals = Marginals .* lognpdf(Scenarios(:, Asset), LogMean(Asset), LogStd(Asset));
end
CopulaDensity = t_copula_density(Uniforms, Correlation, Model.CopulaDF);
Density = CopulaDensity .* Marginals;
Density(Marginals == 0) = 0;            % marginal density underflow (u_i = 0 or 1)
end


function Density = t_copula_density(Uniforms, Correlation, DF)
% t-copula density c(u) = g(x) / (t(x_1) t(x_2)), x_i = T^{-1}(u_i), where g is
% the bivariate t density and t, T the univariate t density and CDF. Same
% values as copulapdf('t', ...), but the quantiles come from betaincinv on the
% smaller tail, so they stay finite for uniforms near 0 or 1 (tinv returns
% -Inf below about 1e-17, and copulapdf then returns NaN).
% formula: T^{-1}(p) = -sqrt(DF (1/z - 1)), z = I^{-1}_{2p}(DF/2, 1/2), p < 1/2
Tail = min(Uniforms, 1 - Uniforms);                     % smaller tail probability
Beta = betaincinv(2 * Tail, DF / 2, 0.5);
Quantiles = sqrt(DF * (1 ./ Beta - 1));
Quantiles(Uniforms < 0.5) = -Quantiles(Uniforms < 0.5);
Joint = mvtpdf(Quantiles, Correlation, DF);             % bivariate t density
Univariate = tpdf(Quantiles(:, 1), DF) .* tpdf(Quantiles(:, 2), DF);
Density = Joint ./ Univariate;
end
