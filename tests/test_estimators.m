function MaxDifference = test_estimators()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test_estimators.m                                                           %
% Checks catbond_one_path.m and credit_one_path.m against a direct            %
% implementation of Equations 2, 6, 9 and 15 that stores every output and     %
% computes every likelihood ratio, on the same random numbers.               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output(s):
% -- MaxDifference: largest absolute difference (error if above 1e-12)

CodeFolder = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(CodeFolder, 'functions'));
MaxDifference = 0;

% ===== CAT bond: 3 macro replications x 12 periods x 50 replications
Model = catbond_settings();
for Macro = 1:3
    Stream = RandStream('mrg32k3a', 'Seed', 21); Stream.Substream = Macro;
    Fast = catbond_one_path(Model, 12, 50, Stream, false);
    Stream = RandStream('mrg32k3a', 'Seed', 21); Stream.Substream = Macro;
    Parameters = catbond_states(Model, 12, Stream);
    Losses = cell(12, 1); Payoffs = cell(12, 1);
    for n = 1:12
        [Losses{n}, Payoffs{n}] = catbond_simulate(Model, Parameters(n, :), 50, Stream);
        Direct = direct_estimates(Losses, Payoffs, Parameters, n, 50, Parameters(n, :), @catbond_density_wrapper);
        MaxDifference = max(MaxDifference, max(abs(Fast.Estimates.current(n, 1:4) - Direct)));
        Direct = direct_estimates(Losses, Payoffs, Parameters, n, 50, Model.FixedStates.hi, @catbond_density_wrapper);
        MaxDifference = max(MaxDifference, max(abs(Fast.Estimates.hi(n, 2:4) - Direct(2:4))));
    end
end

% ===== Credit risk: 2 macro replications x 8 periods x 300 replications
Model = credit_settings();
for Macro = 1:2
    Stream = RandStream('mrg32k3a', 'Seed', 22); Stream.Substream = Macro;
    Fast = credit_one_path(Model, 8, 300, Stream);
    Stream = RandStream('mrg32k3a', 'Seed', 22); Stream.Substream = Macro;
    States = credit_states(Model, 8, Stream);
    Scenarios = cell(8, 1); Outputs = cell(8, 1);
    Density = @(Y, State) credit_density(Model, Y, State);
    for n = 1:8
        Scenarios{n} = credit_scenarios(Model, States(n, :), 300, Stream);
        Outputs{n} = double(credit_loss(Model, Scenarios{n}) > Model.Threshold);
        Direct = direct_estimates(Scenarios, Outputs, States, n, 300, States(n, :), Density);
        MaxDifference = max(MaxDifference, max(abs(Fast.Estimates(n, :) - Direct)));
    end
end
% ===== Underflow: an output with zero target density (and, as happens far
% in the tails, zero mixture density) must contribute 0, not NaN
[I, E, M] = lr_estimators([1; 1; 1], [1; 2; 2], [0.5; 0; 0.2], [0.3; 0; 0.4], [0.6; 0; 0.1], 2, 2);
Expected = [mean([1 * 0.6 / 0.5 / 2, 1 * 0.1 / 0.2 / 2]), NaN, (0.6 / 0.3 + 0.1 / 0.4) / 4];
assert(all(isfinite([I, E, M])), 'LR estimators are not finite when densities underflow');
MaxDifference = max(MaxDifference, max(abs([I, M] - Expected([1 3]))));
fprintf('test_estimators: largest difference from the direct formulas %.2e\n', MaxDifference);
if MaxDifference > 1e-12
    error('The estimators differ from Equations 2, 6, 9 and 15.');
end
end


function Estimates = direct_estimates(Inputs, Outputs, States, n, r, Target, Density)
% SMC, ILR, EWLR and MLR at a target state, with all outputs of experiments
% 1, ..., n (Inputs{i}: random inputs Y of experiment i; Outputs{i}: F(Y))
BatchMeans = zeros(n, 1); BatchVariances = zeros(n, 1); MLRSum = 0;
for i = 1:n
    Y = Inputs{i}; F = Outputs{i};
    Ratio = Density(Y, Target) ./ Density(Y, States(i, :));          % h(Y; x)/h(Y; x_i)
    Mixture = zeros(size(F));
    for l = 1:n
        Mixture = Mixture + Density(Y, States(l, :)) / n;              % h*_n(Y)
    end
    BatchMeans(i) = mean(F .* Ratio);
    BatchVariances(i) = var(F .* Ratio);
    MLRSum = MLRSum + sum(F .* Density(Y, Target) ./ Mixture);
end
InverseVariances = zeros(n, 1);
InverseVariances(BatchVariances > 0) = 1 ./ BatchVariances(BatchVariances > 0);
if sum(InverseVariances) > 0
    EWLR = sum(InverseVariances .* BatchMeans) / sum(InverseVariances);
else
    EWLR = 0;
end
Estimates = [mean(Outputs{n}), mean(BatchMeans), EWLR, MLRSum / (n * r)];
end


function Density = catbond_density_wrapper(Losses, Parameters)
% catbond_density with the argument order of direct_estimates
Density = catbond_density(Losses, Parameters);
end
