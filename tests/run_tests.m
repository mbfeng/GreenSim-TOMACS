function run_tests()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run_tests.m                                                                 %
% Regression tests: (1) the estimators equal a direct implementation of the   %
% paper's formulas; (2) stochastic kriging and credit density checks;         %
% (3) all outputs equal the saved snapshot.                                   %
% Run from the Code folder: addpath('tests'); run_tests                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TestsFolder = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(TestsFolder), 'functions'), TestsFolder);

test_estimators();
test_sk();
test_credit_density();

Saved = load(fullfile(TestsFolder, 'reference_outputs.mat'));
Saved = Saved.Reference;
Current = reference_outputs();
Differences = [max(abs(Current.Prices - Saved.Prices));
               max(abs(Current.Probabilities - Saved.Probabilities));
               max(abs(Current.Published - Saved.Published))];
for Macro = 1:2
    Both = ~isnan(Saved.CatBond{Macro});
    Differences(end + 1) = max(abs(Current.CatBond{Macro}(Both) - Saved.CatBond{Macro}(Both))); %#ok<AGROW>
    Differences(end + 1) = max(abs(Current.Credit{Macro} - Saved.Credit{Macro}), [], 'all'); %#ok<AGROW>
end
fprintf('run_tests: largest difference from the snapshot %.2e\n', max(Differences));
if max(Differences) > 1e-8
    error('The results differ from tests/reference_outputs.mat.');
end
fprintf('run_tests: all tests passed\n');
end
