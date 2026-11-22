function extract_published_runs(CatBondFile, CreditFile, RunFolder)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract_published_runs.m                                                    %
% One-off provenance tool: reduces the saved workspaces of the runs behind    %
% Figures 3-5 (CAT bond: results_10KReps_with_Trend_EB.mat, 2016; credit      %
% risk: results_with_GIBF.mat, 2016) to per-period summaries                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input(s):
% -- CatBondFile, CreditFile: paths of the saved workspaces (not in this
%    repository)
% -- RunFolder: output folder, e.g. results/published_run
% Output(s):
% -- RunFolder/catbond_summary.csv: Period, Target (current, lo, mi, hi),
%    Estimator (SMC, ILR, EWLR, MLR, GSK), MSE, SE, Top1Share
% -- RunFolder/credit_summary.csv: Period, Target (current), Estimator (SMC,
%    ILR, EWLR, MLR), MSE, SE, Top1Share
% MSE = mean over the 10,000 macro replications of the squared error, SE =
% its standard error (std/sqrt(10,000)), Top1Share = share of the MSE from
% the worst 1% of macro replications (NaN where only MSE and SE were saved).
% Names in the original code: WLR = EWLR, SK1 = GSK. The fixed-state SMC
% variances were computed exactly (r = 100 replications), so their SE is 0.

if ~exist(RunFolder, 'dir'), mkdir(RunFolder); end
Estimators = {'SMC', 'ILR', 'EWLR', 'MLR', 'GSK'};
CodeNames = {'SMC', 'ILR', 'WLR', 'MLR', 'SK1'};

% ===== CAT bond
Saved = load(CatBondFile);
Rows = {};
for e = 1:numel(Estimators)
    % current state: per-macro squared errors are saved
    Errors = Saved.(['RMSE_' CodeNames{e}]);            % 100 x 10,000
    Rows = [Rows; summary_rows('current', Estimators{e}, Errors)]; %#ok<AGROW>
    % fixed states: only mean and standard error were saved
    for Target = {'lo', 'mi', 'hi'}
        Name = Target{1};
        MSE = Saved.(['RMSE_Y_' CodeNames{e} '_' Name]);
        if strcmp(Estimators{e}, 'SMC')
            SE = zeros(size(MSE));                       % exact variance
        else
            SE = Saved.(['RMSE_E_' CodeNames{e} '_' Name]);
        end
        Rows = [Rows; mean_rows(Name, Estimators{e}, MSE, SE)]; %#ok<AGROW>
    end
end
write_summary(Rows, fullfile(RunFolder, 'catbond_summary.csv'));
clear Saved

% ===== Credit risk
Saved = load(CreditFile, 'RMSE_SMC', 'RMSE_ILR', 'RMSE_WLR', 'RMSE_MLR');
Rows = {};
for e = 1:4
    Errors = Saved.(['RMSE_' CodeNames{e}]);            % 52 x 10,000
    Rows = [Rows; summary_rows('current', Estimators{e}, Errors)]; %#ok<AGROW>
end
write_summary(Rows, fullfile(RunFolder, 'credit_summary.csv'));
fprintf('extract_published_runs: summaries written to %s\n', RunFolder);
end


function Rows = summary_rows(Target, Estimator, Errors)
% One row per period from per-macro squared errors (periods x macros)
[NumPeriods, NumMacro] = size(Errors);
Rows = cell(NumPeriods, 6);
for n = 1:NumPeriods
    Sorted = sort(Errors(n, :), 'descend');
    NumTop = max(1, round(0.01 * NumMacro));
    Rows(n, :) = {n, Target, Estimator, mean(Errors(n, :)), ...
        std(Errors(n, :)) / sqrt(NumMacro), sum(Sorted(1:NumTop)) / sum(Sorted)};
end
end


function Rows = mean_rows(Target, Estimator, MSE, SE)
% One row per period from saved means and standard errors
NumPeriods = numel(MSE);
Rows = cell(NumPeriods, 6);
for n = 1:NumPeriods
    Rows(n, :) = {n, Target, Estimator, MSE(n), SE(n), NaN};
end
end


function write_summary(Rows, File)
% Writes the summary rows as CSV (17 significant digits)
Table = cell2table(Rows, 'VariableNames', {'Period', 'Target', 'Estimator', 'MSE', 'SE', 'Top1Share'});
% GSK is fitted from period 5 on; earlier rows hold no estimates
Table = Table(~(strcmp(Table.Estimator, 'GSK') & Table.Period < 5), :);
writetable(Table, File);
end
