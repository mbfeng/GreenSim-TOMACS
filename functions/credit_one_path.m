function Result = credit_one_path(Model, NumPeriods, NumReplications, Stream)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% credit_one_path.m                                                           %
% One macro replication of the credit risk example (Section 6.2): a weekly   %
% sample path of asset values, the n-th experiment's r replications, and the %
% SMC, ILR, EWLR and MLR estimates at the current state                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input(s):
% -- Model: struct from credit_settings.m
% -- NumPeriods: N (52 in the paper); NumReplications: r (1,000)
% -- Stream: RandStream of this macro replication
% Output(s): struct Result with
% -- .Estimates: N x 4, columns SMC, ILR, EWLR, MLR
% -- .TrueProbabilities: N x 1, mu(x_n)
% Outputs are F(Y) = 1{L_1 + L_2 > kappa}; only outputs with F = 1 are stored,
% since the others add nothing to the LR estimators' sums.

States = credit_states(Model, NumPeriods, Stream);             % x_1 = (100, 90)
TrueProbabilities = credit_true_probability(Model, States);
Estimates = zeros(NumPeriods, 4);

Stored = zeros(0, 2);                        % scenarios with F = 1
Batch = zeros(0, 1);                         % experiment of each
Density = zeros(0, NumPeriods);              % Density(j, i) = h(Y_j; x_i)
for n = 1:NumPeriods
    Scenarios = credit_scenarios(Model, States(n, :), NumReplications, Stream);
    HighLoss = credit_loss(Model, Scenarios) > Model.Threshold;   % F(Y) = 1
    Estimates(n, 1) = mean(HighLoss);                              % SMC

    NewScenarios = Scenarios(HighLoss, :);
    NumNew = size(NewScenarios, 1);
    NewRows = size(Stored, 1) + (1:NumNew);
    Stored = [Stored; NewScenarios];                               %#ok<AGROW>
    Batch = [Batch; n * ones(NumNew, 1)];                          %#ok<AGROW>
    Density = [Density; zeros(NumNew, NumPeriods)];                %#ok<AGROW>
    for i = 1:n - 1
        Density(NewRows, i) = credit_density(Model, NewScenarios, States(i, :));
    end
    Density(:, n) = credit_density(Model, Stored, States(n, :));

    OwnDensity = Density(sub2ind(size(Density), (1:size(Stored, 1))', Batch));
    MixtureDensity = mean(Density(:, 1:n), 2);
    Outputs = ones(size(Stored, 1), 1);
    [Estimates(n, 2), Estimates(n, 3), Estimates(n, 4)] = lr_estimators(Outputs, Batch, ...
        OwnDensity, MixtureDensity, Density(:, n), n, NumReplications);
end
Result.Estimates = Estimates;
Result.TrueProbabilities = TrueProbabilities;
end
