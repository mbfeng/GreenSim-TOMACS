function Result = catbond_one_path(Model, NumPeriods, NumReplications, Stream, UseGSK)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% catbond_one_path.m                                                          %
% One macro replication of the CAT bond example (Section 6.1): a sample path  %
% of states, the n-th experiment's r replications, and the SMC, ILR, EWLR,    %
% MLR and GSK estimates at the current state and at three fixed states,       %
% all from the same simulation output (common random numbers)                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input(s):
% -- Model: struct from catbond_settings.m
% -- NumPeriods: N (100 in the paper); NumReplications: r (100)
% -- Stream: RandStream of this macro replication
% -- UseGSK: true to compute the green stochastic kriging estimator
% Output(s): struct Result with
% -- .Estimates: struct with fields current, lo, mi, hi, each N x 5,
%    columns SMC, ILR, EWLR, MLR, GSK (GSK is NaN before period 5; the
%    fixed-state SMC column is NaN: its variance is known exactly)
% -- .TruePrices: N x 1, mu(x_n) of the current states
% Algorithms 1 (ILR, EWLR) and 2 (MLR): densities of every stored output
% under every state so far are kept in a matrix filled period by period.

GSKStart = 5;                                     % first period with a GSK fit
Parameters = catbond_states(Model, NumPeriods, Stream);     % (lambda_n, theta_n)
TruePrices = catbond_price(Model, Parameters);
FixedNames = {'lo', 'mi', 'hi'};

NumOutputs = NumPeriods * NumReplications;
Losses = zeros(NumOutputs, 1);                    % Y of all experiments
Payoffs = zeros(NumOutputs, 1);                   % F(Y)
Batch = zeros(NumOutputs, 1);                     % experiment of each output
Density = zeros(NumOutputs, NumPeriods);          % Density(j, i) = h(Y_j; x_i)
FixedDensity = zeros(NumOutputs, 3);              % h(Y_j; x_lo), h(Y_j; x_mi), h(Y_j; x_hi)
SMCMean = zeros(NumPeriods, 1);                   % responses for GSK
SMCVariance = zeros(NumPeriods, 1);               % intrinsic variances for GSK

Estimates.current = nan(NumPeriods, 5);
for f = 1:3
    Estimates.(FixedNames{f}) = nan(NumPeriods, 5);
end
FitGSK = [];

for n = 1:NumPeriods
    %-- Experiment n: r replications at the current parameters
    New = (n - 1) * NumReplications + (1:NumReplications);      % rows of the new outputs
    [Losses(New), Payoffs(New)] = catbond_simulate(Model, Parameters(n, :), NumReplications, Stream);
    Batch(New) = n;
    Stored = 1:n * NumReplications;                              % outputs so far

    %-- Densities: new outputs under the earlier states, all outputs under x_n
    for i = 1:n - 1
        Density(New, i) = catbond_density(Losses(New), Parameters(i, :));
    end
    Density(Stored, n) = catbond_density(Losses(Stored), Parameters(n, :));
    for f = 1:3
        FixedDensity(New, f) = catbond_density(Losses(New), Model.FixedStates.(FixedNames{f}));
    end
    OwnDensity = Density(sub2ind(size(Density), Stored', Batch(Stored)));
    MixtureDensity = mean(Density(Stored, 1:n), 2);

    %-- SMC at the current state
    SMCMean(n) = mean(Payoffs(New));
    SMCVariance(n) = var(Payoffs(New)) / NumReplications;
    Estimates.current(n, 1) = SMCMean(n);

    %-- LR estimators at the current state and at the fixed states
    [Estimates.current(n, 2), Estimates.current(n, 3), Estimates.current(n, 4)] = lr_estimators( ...
        Payoffs(Stored), Batch(Stored), OwnDensity, MixtureDensity, Density(Stored, n), n, NumReplications);
    for f = 1:3
        [I, E, M] = lr_estimators(Payoffs(Stored), Batch(Stored), OwnDensity, MixtureDensity, ...
            FixedDensity(Stored, f), n, NumReplications);
        Estimates.(FixedNames{f})(n, 2:4) = [I, E, M];
    end

    %-- Green stochastic kriging: SMC estimates of all experiments so far
    % as responses, warm-started from the previous period's fit
    if UseGSK && n >= GSKStart
        FitGSK = sk_fit(Parameters(1:n, :), SMCMean(1:n), SMCVariance(1:n), FitGSK);
        Estimates.current(n, 5) = sk_predict(FitGSK, Parameters(n, :));
        for f = 1:3
            Estimates.(FixedNames{f})(n, 5) = sk_predict(FitGSK, Model.FixedStates.(FixedNames{f}));
        end
    end
end
Result.Estimates = Estimates;
Result.TruePrices = TruePrices;
end
