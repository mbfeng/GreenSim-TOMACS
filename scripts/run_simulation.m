function run_simulation(Example, NumMacro, Seed, NumWorkers, RunFolder, NumParts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run_simulation.m                                                            %
% Macro replications of an example of Section 6 and the summary of the MSEs  %
% of its estimators (Figures 3-5)                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input(s) (all but Example optional; defaults are the paper's settings):
% -- Example: 'catbond' (Section 6.1: 100 periods, r = 100, with GSK) or
%    'credit' (Section 6.2: 52 periods, r = 1,000)
% -- NumMacro: macro replications, 10,000
% -- Seed: seed of the random streams, 2017
% -- NumWorkers: parallel workers (0 = serial), 0
% -- RunFolder: output folder, results/rerun_seed<Seed>
% -- NumParts: macro replications are run and saved in parts, so an
%    interrupted run can be resumed (finished parts are skipped), 20
% Output(s): files in RunFolder
% -- parts/<Example>_part_<p>.mat: per-macro squared errors of the part
% -- <Example>_summary.csv: Period, Target, Estimator, MSE, SE, Top1Share
%    (as in results/published_run); run_info_<Example>.txt
% Reproducibility: macro replication m uses substream m of one mrg32k3a
% stream with the given seed; serial and parallel runs are identical.
% Example (server, from the Code folder):
%   matlab -batch "addpath('functions','scripts'); run_simulation('catbond', 10000, 2017, 40)"

Example = char(Example);                      % accept strings as well as chars
if nargin < 2 || isempty(NumMacro), NumMacro = 10000; end
if nargin < 3 || isempty(Seed), Seed = 2017; end
if nargin < 4 || isempty(NumWorkers), NumWorkers = 0; end
if nargin < 5 || isempty(RunFolder), RunFolder = fullfile('results', sprintf('rerun_seed%d', Seed)); end
if nargin < 6 || isempty(NumParts), NumParts = 20; end
RunFolder = char(RunFolder);
switch Example
    case 'catbond'
        Model = catbond_settings(); NumPeriods = 100; NumReplications = 100;
        Targets = {'current', 'lo', 'mi', 'hi'};
        Estimators = {'SMC', 'ILR', 'EWLR', 'MLR', 'GSK'};
    case 'credit'
        Model = credit_settings(); NumPeriods = 52; NumReplications = 1000;
        Targets = {'current'};
        Estimators = {'SMC', 'ILR', 'EWLR', 'MLR'};
    otherwise
        error('Example must be ''catbond'' or ''credit''.');
end
PartsFolder = fullfile(RunFolder, 'parts');
if ~exist(PartsFolder, 'dir'), mkdir(PartsFolder); end
if NumWorkers > 0 && isempty(gcp('nocreate'))
    Cluster = parcluster('Processes');              % allow the requested number
    Cluster.NumWorkers = max(Cluster.NumWorkers, NumWorkers);   % of workers
    parpool(Cluster, NumWorkers);
end

% ===== Parts: squared errors, periods x estimators x targets x macros
Start = tic;
PerPart = ceil(NumMacro / NumParts);
for Part = 1:NumParts
    PartFile = fullfile(PartsFolder, sprintf('%s_part_%03d.mat', Example, Part));
    if exist(PartFile, 'file'), continue; end
    Macros = (Part - 1) * PerPart + 1 : min(Part * PerPart, NumMacro);
    if isempty(Macros), break; end
    Errors = zeros(NumPeriods, numel(Estimators), numel(Targets), numel(Macros), 'single');
    parfor (m = 1:numel(Macros), max(NumWorkers, 0))
        Stream = RandStream('mrg32k3a', 'Seed', Seed);
        Stream.Substream = Macros(m);
        Errors(:, :, :, m) = one_macro(Example, Model, NumPeriods, NumReplications, Stream, Targets);
    end
    save(PartFile, 'Errors', 'Macros', '-v7.3');
    fprintf('%s: part %d of %d done, %.1f min so far\n', Example, Part, NumParts, toc(Start) / 60);
end

% ===== Summary over all parts
All = cell(NumParts, 1);
for Part = 1:NumParts
    PartFile = fullfile(PartsFolder, sprintf('%s_part_%03d.mat', Example, Part));
    if exist(PartFile, 'file')
        Loaded = load(PartFile, 'Errors');
        All{Part} = double(Loaded.Errors);
    end
end
Errors = cat(4, All{:});                                   % periods x est x targets x macros
NumDone = size(Errors, 4);
Rows = {};
for t = 1:numel(Targets)
    for e = 1:numel(Estimators)
        for n = 1:NumPeriods
            if strcmp(Example, 'catbond') && ~strcmp(Targets{t}, 'current') && strcmp(Estimators{e}, 'SMC')
                [~, Variance] = catbond_price(Model, Model.FixedStates.(Targets{t}), NumReplications);
                Rows(end + 1, :) = {n, Targets{t}, Estimators{e}, Variance, 0, NaN}; %#ok<AGROW>
                continue                                   % exact SMC variance
            end
            Values = squeeze(Errors(n, e, t, :));
            if all(isnan(Values)), continue; end           % GSK before period 5
            Sorted = sort(Values, 'descend');
            NumTop = max(1, round(0.01 * NumDone));
            Rows(end + 1, :) = {n, Targets{t}, Estimators{e}, mean(Values), ...
                std(Values) / sqrt(NumDone), sum(Sorted(1:NumTop)) / sum(Sorted)}; %#ok<AGROW>
        end
    end
end
Summary = cell2table(Rows, 'VariableNames', {'Period', 'Target', 'Estimator', 'MSE', 'SE', 'Top1Share'});
writetable(Summary, fullfile(RunFolder, [Example '_summary.csv']));
Info = fopen(fullfile(RunFolder, ['run_info_' Example '.txt']), 'w');
fprintf(Info, '%s: macro replications %d, periods %d, replications %d, seed %d\n', ...
    Example, NumDone, NumPeriods, NumReplications, Seed);
fprintf(Info, 'workers %d, wall-clock of this call %.1f min, MATLAB %s\n', NumWorkers, toc(Start) / 60, version);
fclose(Info);
end


function Errors = one_macro(Example, Model, NumPeriods, NumReplications, Stream, Targets)
% Squared errors of all estimators in one macro replication:
% periods x estimators x targets
if strcmp(Example, 'catbond')
    Result = catbond_one_path(Model, NumPeriods, NumReplications, Stream, true);
    Errors = zeros(NumPeriods, 5, numel(Targets));
    Errors(:, :, 1) = (Result.Estimates.current - Result.TruePrices).^2;
    for t = 2:numel(Targets)
        TruePrice = catbond_price(Model, Model.FixedStates.(Targets{t}));
        Errors(:, :, t) = (Result.Estimates.(Targets{t}) - TruePrice).^2;
    end
else
    Result = credit_one_path(Model, NumPeriods, NumReplications, Stream);
    Errors = (Result.Estimates - Result.TrueProbabilities).^2;
end
end
