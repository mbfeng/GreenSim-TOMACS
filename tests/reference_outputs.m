function Reference = reference_outputs()
% Outputs on fixed inputs and seeds, used by make_reference_outputs.m (saves
% them) and run_tests.m (compares with them)
% Output(s): struct Reference with
% -- .CatBond: estimates of 2 macro replications (30 periods, r = 100, GSK)
% -- .Credit: estimates of 2 macro replications (10 periods, r = 500)
% -- .Prices, .Probabilities: true values at fixed states
% -- .Published: summaries of the published runs (as read by the figures)
CatBondModel = catbond_settings();
CreditModel = credit_settings();
Reference.CatBond = cell(2, 1);
Reference.Credit = cell(2, 1);
for Macro = 1:2
    Stream = RandStream('mrg32k3a', 'Seed', 51); Stream.Substream = Macro;
    Result = catbond_one_path(CatBondModel, 30, 100, Stream, true);
    Reference.CatBond{Macro} = [Result.Estimates.current, Result.Estimates.mi, Result.Estimates.hi];
    Stream = RandStream('mrg32k3a', 'Seed', 52); Stream.Substream = Macro;
    Result = credit_one_path(CreditModel, 10, 500, Stream);
    Reference.Credit{Macro} = Result.Estimates;
end
Reference.Prices = catbond_price(CatBondModel, [2 4; 3 5; 4 6; 2.5 5.5]);
Reference.Probabilities = credit_true_probability(CreditModel, [100 90; 90 85; 110 95]);
CodeFolder = fileparts(fileparts(mfilename('fullpath')));
Published = fullfile(CodeFolder, 'results', 'published_run');
Reference.Published = [readtable(fullfile(Published, 'catbond_summary.csv')).MSE; ...
                       readtable(fullfile(Published, 'credit_summary.csv')).MSE];
end
