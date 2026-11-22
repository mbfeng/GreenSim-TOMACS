function Report = verify_against_paper(RunFolder, ReportFile)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% verify_against_paper.m                                                      %
% Checks the MSE curves of Figures 3-5 against the published paper            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input(s):
% -- RunFolder: results/published_run (exact mode) or a rerun folder (Monte
%    Carlo mode); each holds catbond_summary.csv and credit_summary.csv
% -- ReportFile: CSV report (none written if empty)
% Output(s):
% -- Report: table, one row per curve (exact mode) or per period, target and
%    estimator (Monte Carlo mode)
% Exact mode: the curves plotted in the paper's vector figures
%   (verify/paper_*_points.csv, read by tools/digitize_published_figures.py)
%   are compared with the published-run summary. The axis tick labels of the
%   figures are drawn as outlines, so for each panel the log-log axis scales
%   (x_px = a + b log10(n), y_px = c + d log10(MSE)) are fitted to all its
%   curves jointly; every plotted point must then lie within 0.002 px of the
%   summary (the coordinates are printed to 3 decimals).
% Monte Carlo mode: z = (x - v) / sqrt(se_x^2 + se_v^2) per period, target
%   and estimator; |z| > 2.5 is flagged. When the worst 1% of macro
%   replications contribute more than half of an MSE (Top1Share > 0.5), its
%   standard error is unreliable and the status is 'heavy-tailed'.

CodeFolder = fileparts(fileparts(mfilename('fullpath')));
IsPublished = endsWith(strip(RunFolder, 'right', filesep), 'published_run');
CatBond = readtable(fullfile(RunFolder, 'catbond_summary.csv'), 'TextType', 'string');
Credit = readtable(fullfile(RunFolder, 'credit_summary.csv'), 'TextType', 'string');

if IsPublished
    Report = exact_mode(CodeFolder, CatBond, Credit);
    NumOk = sum(Report.Status == "ok");
    fprintf('verify_against_paper (exact, published run): %d of %d plotted curves of Figures 3-5 reproduced (largest residual %.4f px).\n', ...
        NumOk, height(Report), max(Report.MaxResidualPx));
else
    Published = fullfile(CodeFolder, 'results', 'published_run');
    Report = [monte_carlo_mode(CatBond, readtable(fullfile(Published, 'catbond_summary.csv'), 'TextType', 'string'), "CAT bond");
              monte_carlo_mode(Credit, readtable(fullfile(Published, 'credit_summary.csv'), 'TextType', 'string'), "credit risk")];
    NumOk = sum(Report.Status == "ok" | Report.Status == "heavy-tailed");
    fprintf('verify_against_paper (Monte Carlo, %s vs published run): %d of %d ok (of which %d heavy-tailed); expected chance flags at |z| > 2.5: %.1f.\n', ...
        RunFolder, NumOk, height(Report), sum(Report.Status == "heavy-tailed"), 0.0124 * height(Report));
end
if ~isempty(ReportFile)
    writetable(Report, ReportFile);
end
Flagged = Report(~ismember(Report.Status, ["ok", "heavy-tailed"]), :);
if height(Flagged) > 0
    disp(Flagged);
end
if IsPublished && NumOk < height(Report)
    error('Some plotted curves were not reproduced.');
end
end


function Report = exact_mode(CodeFolder, CatBond, Credit)
% Exact mode: fit axis scales per panel and compare all plotted curves
Colours = ["0 0.8 0.8", "1 0.4 0.4", "1 0.8 0", "0.8 0.2 1", "0 0 0.8"];
Estimators = ["SMC", "ILR", "EWLR", "MLR", "GSK"];
% {figure file, summary, number of panels, panel targets}
Figures = {"sig2_lomihi", CatBond, ["mi", "hi"], 5;
           "sig2n", CatBond, "current", 5;
           "sig2n_example2", Credit, "current", 4};
Rows = {};
for f = 1:size(Figures, 1)
    Points = readtable(fullfile(CodeFolder, 'verify', "paper_" + Figures{f, 1} + "_points.csv"), 'TextType', 'string');
    Summary = Figures{f, 2};
    Targets = Figures{f, 3};
    NumEstimators = Figures{f, 4};
    for p = 1:numel(Targets)
        % collect (n, MSE, x_px, y_px) for all curves of this panel
        Curves = cell(NumEstimators, 1);
        for e = 1:NumEstimators
            Lines = curve_lines(Points, Colours(e), e == 5);
            Line = Lines{p};                                  % panel p
            Rows_e = Summary(Summary.Target == Targets(p) & Summary.Estimator == Estimators(e), :);
            Periods = Rows_e.Period;
            if numel(Periods) ~= height(Line)
                error('Curve %s/%s/%s: %d plotted points, %d periods.', Figures{f, 1}, Targets(p), Estimators(e), height(Line), numel(Periods));
            end
            Curves{e} = [Periods, Rows_e.MSE, Line.X, Line.Y];
        end
        All = vertcat(Curves{:});
        X = [ones(size(All, 1), 1), log10(All(:, 1))];
        XFit = X \ All(:, 3);                                % x_px = a + b log10(n)
        Y = [ones(size(All, 1), 1), log10(All(:, 2))];
        YFit = Y \ All(:, 4);                                % y_px = c + d log10(MSE)
        for e = 1:NumEstimators
            C = Curves{e};
            Residual = max([abs([ones(size(C, 1), 1), log10(C(:, 1))] * XFit - C(:, 3));
                            abs([ones(size(C, 1), 1), log10(C(:, 2))] * YFit - C(:, 4))]);
            if Residual <= 0.002
                Status = "ok";
            else
                Status = "UNEXPECTED";
            end
            Rows(end + 1, :) = {Figures{f, 1}, Targets(p), Estimators(e), size(C, 1), Residual, Status}; %#ok<AGROW>
        end
    end
end
Report = cell2table(Rows, 'VariableNames', {'Figure', 'Target', 'Estimator', 'NumPoints', 'MaxResidualPx', 'Status'});
end


function Lines = curve_lines(Points, Colour, IsDashedGSK)
% Polylines of one colour that are curves (not legend samples, error bars or
% the slope -1 reference line), in the order drawn (= panel order)
Candidates = unique(Points.Line(Points.Colour == Colour), 'stable');
Lines = {};
for k = 1:numel(Candidates)
    Line = Points(Points.Line == Candidates(k), :);
    IsCurve = height(Line) >= 40;                     % curves have 52-100 points
    if IsDashedGSK
        IsCurve = IsCurve && Line.Dash(1) ~= "";      % GSK is dashed; the
    end                                               % reference line is solid
    if IsCurve
        Lines{end + 1} = Line; %#ok<AGROW>
    end
end
end


function Report = monte_carlo_mode(Run, Published, Example)
% Monte Carlo mode: z-test per period, target and estimator
Rows = {};
for i = 1:height(Run)
    Match = Published(Published.Period == Run.Period(i) & Published.Target == Run.Target(i) ...
        & Published.Estimator == Run.Estimator(i), :);
    if height(Match) ~= 1, continue; end
    Denominator = sqrt(Run.SE(i)^2 + Match.SE^2);
    if Denominator == 0
        Z = (Run.MSE(i) - Match.MSE) / max(abs(Match.MSE), eps) * 1e6;   % exact values
    else
        Z = (Run.MSE(i) - Match.MSE) / Denominator;
    end
    HeavyTailed = max(Run.Top1Share(i), Match.Top1Share) > 0.5;
    if abs(Z) <= 2.5
        Status = "ok";
    elseif HeavyTailed
        Status = "heavy-tailed";
    else
        Status = "UNEXPECTED";
    end
    Rows(end + 1, :) = {Example, Run.Period(i), Run.Target(i), Run.Estimator(i), ...
        Run.MSE(i), Match.MSE, Z, Run.Top1Share(i), Status}; %#ok<AGROW>
end
Report = cell2table(Rows, 'VariableNames', {'Example', 'Period', 'Target', 'Estimator', ...
    'Rerun', 'Published', 'Z', 'Top1Share', 'Status'});
end
