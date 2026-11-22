%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reproduce_paper.m                                                           %
% Entry point: Figures 2-5 of Feng and Staum (2017), "Green simulation:       %
% Reusing the output of repeated experiments", ACM Transactions on Modeling  %
% and Computer Simulation 27(4), Article 23, and the check of Figures 3-5     %
% against the published paper. Run from this folder.                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ===== Settings
RERUN_SIMULATIONS = false;  % false: draw the figures from the saved summaries
                            % of the runs behind the published figures
                            % (seconds); true: rerun the 2 x 10,000 macro
                            % replications first (hours; see README)
NUM_WORKERS = 0;            % parallel workers for a rerun (0 = serial)
SEED = 2017;                % seed of a rerun

% ===== Paths
CodeFolder = fileparts(mfilename('fullpath'));
addpath(fullfile(CodeFolder, 'functions'), fullfile(CodeFolder, 'scripts'), ...
    fullfile(CodeFolder, 'verify'));
cd(CodeFolder);

% ===== Runs
if RERUN_SIMULATIONS
    RunFolder = fullfile('results', sprintf('rerun_seed%d', SEED));
    run_simulation('catbond', 10000, SEED, NUM_WORKERS, RunFolder);
    run_simulation('credit', 10000, SEED, NUM_WORKERS, RunFolder);
else
    RunFolder = fullfile('results', 'published_run');
end

% ===== Figures 2-5 (Figure 1 is drawn in the LaTeX source)
make_figures(RunFolder, 'figures');
fprintf('Figures 2-5 written to figures/ (from %s)\n', RunFolder);

% ===== Check against the published figures
verify_against_paper(RunFolder, fullfile('verify', ['verification_' strrep(RunFolder, filesep, '_') '.csv']));
