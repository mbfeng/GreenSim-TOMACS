function MaxDifference = test_credit_density()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test_credit_density.m                                                       %
% credit_density.m equals the t-copula density of copulapdf times the        %
% lognormal marginals in the bulk, and stays finite and positive for extreme %
% scenarios, where tinv (inside copulapdf) returns -Inf and copulapdf NaN    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Model = credit_settings();
State = [100, 90];
Stream = RandStream('mrg32k3a', 'Seed', 1);
Scenarios = 100 * exp(0.1 * randn(Stream, 2000, 2));       % bulk of the distribution
[LogMean, LogStd] = scenario_marginals(Model, State);
Uniforms = [logncdf(Scenarios(:, 1), LogMean(1), LogStd(1)), logncdf(Scenarios(:, 2), LogMean(2), LogStd(2))];
Correlation = [1, Model.CopulaRho; Model.CopulaRho, 1];
Direct = copulapdf('t', Uniforms, Correlation, Model.CopulaDF) .* ...
    lognpdf(Scenarios(:, 1), LogMean(1), LogStd(1)) .* lognpdf(Scenarios(:, 2), LogMean(2), LogStd(2));
MaxDifference = max(abs(credit_density(Model, Scenarios, State) - Direct) ./ Direct);

Extreme = exp(LogMean - 9 * LogStd);                        % u_i of order 1e-19
ExtremeDensity = credit_density(Model, Extreme, State);
fprintf('test_credit_density: largest relative difference from copulapdf %.2e; extreme scenario density %.3g\n', ...
    MaxDifference, ExtremeDensity);
if MaxDifference > 1e-8 || ~isfinite(ExtremeDensity) || ExtremeDensity <= 0
    error('credit_density.m differs from copulapdf or is not finite in the tails.');
end
end
