function [ILR, EWLR, MLR] = lr_estimators(Outputs, Batch, OwnDensity, MixtureDensity, TargetDensity, NumPeriods, NumReplications)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lr_estimators.m                                                             %
% The three likelihood ratio (LR) estimators of Section 2 for one target      %
% state, from the stored outputs of experiments 1, ..., n                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input(s) (one entry per stored output; outputs with F = 0 may be omitted,
% because they add nothing to the sums below):
% -- Outputs: F(Y_i^(j))
% -- Batch: experiment i in which each output was generated
% -- OwnDensity: h(Y_i^(j); x_i), density under the state that generated it
% -- MixtureDensity: (1/n) sum_{l<=n} h(Y_i^(j); x_l)
% -- TargetDensity: h(Y_i^(j); x), x = the target state
% -- NumPeriods: n; NumReplications: r (outputs per experiment)
% Output(s):
% -- ILR: (1/n) sum_i (1/r) sum_j F h(.; x)/h(.; x_i)                (Eq. 6)
% -- EWLR: sum_i w_i (1/r) sum_j F h(.; x)/h(.; x_i), w_i proportional to
%    the inverse sample variance of experiment i's LR terms           (Eq. 9)
% -- MLR: (1/(nr)) sum_i sum_j F h(.; x) / ((1/n) sum_l h(.; x_l))   (Eq. 15)
% The per-experiment means and sample variances are computed from sums and
% sums of squares, so they are those of all r terms (zeros included).

r = NumReplications;
% An output whose density under the target state is 0 contributes 0 to all
% estimators; computing 0/0 would give NaN when the other densities also
% underflow to 0 (scenarios far in the tails), so such terms are skipped
Contributes = TargetDensity > 0;
Outputs = Outputs(Contributes); Batch = Batch(Contributes);
OwnDensity = OwnDensity(Contributes); MixtureDensity = MixtureDensity(Contributes);
TargetDensity = TargetDensity(Contributes);
Terms = Outputs .* TargetDensity ./ OwnDensity;          % F h(.; x)/h(.; x_i)
Sums = accumarray(Batch(:), Terms(:), [NumPeriods, 1]);           % per experiment
SumSquares = accumarray(Batch(:), Terms(:).^2, [NumPeriods, 1]);
BatchMeans = Sums / r;                                            % ILR of experiment i
BatchVariances = (SumSquares - Sums.^2 / r) / (r - 1);            % sample variances

ILR = mean(BatchMeans);

% EWLR weights: inverse variances, normalized. As in the code that produced
% the published credit risk results, an experiment whose LR terms all equal
% (e.g., no high losses: all terms 0) has zero sample variance and gets
% weight 0; if all experiments have zero variance, the estimate is 0.
InverseVariances = zeros(NumPeriods, 1);
Positive = BatchVariances > 0;
InverseVariances(Positive) = 1 ./ BatchVariances(Positive);
if sum(InverseVariances) > 0
    Weights = InverseVariances / sum(InverseVariances);
else
    Weights = zeros(NumPeriods, 1);
end
EWLR = sum(Weights .* BatchMeans);

MLR = sum(Outputs .* TargetDensity ./ MixtureDensity) / (NumPeriods * r);
end
