function make_figures(RunFolder, FigureFolder)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make_figures.m                                                              %
% Figures 2-5 of the paper, with the file names of its LaTeX source:          %
% paramhist.pdf (Fig. 2), sig2_lomihi.pdf (Fig. 3), sig2n.pdf (Fig. 4),      %
% sig2n_example2.pdf (Fig. 5). Figure 1 is drawn in the LaTeX source.        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input(s):
% -- RunFolder: folder with catbond_summary.csv and credit_summary.csv
%    (results/published_run or a rerun)
% -- FigureFolder: output folder

if ~exist(FigureFolder, 'dir'), mkdir(FigureFolder); end
CatBond = readtable(fullfile(RunFolder, 'catbond_summary.csv'), 'TextType', 'string');
Credit = readtable(fullfile(RunFolder, 'credit_summary.csv'), 'TextType', 'string');

Style.Estimators = ["SMC", "ILR", "EWLR", "MLR", "GSK"];
Style.Lines = {'-', '-.', '-.', '-.', '--'};
Style.Colours = [0 204 204; 255 102 102; 255 204 0; 204 51 255; 0 0 204] / 255;

% ===== Figure 2: parameters sampled from the stationary distribution
Model = catbond_settings();
Stream = RandStream('mrg32k3a', 'Seed', 2);
StationaryStd = Model.NoiseStd ./ sqrt(1 - Model.AutoRegression.^2);
States = (Model.LongRunMean ./ (1 - Model.AutoRegression)) + StationaryStd .* randn(Stream, 1e6, 2);
Parameters = catbond_parameters(Model, States);
Figure = figure('Visible', 'off', 'Position', [100 100 900 420]);
Titles = {'Normalized histogram of $\lambda$', 'Normalized histogram of $\theta$'};
Limits = {Model.LambdaRange, Model.ThetaRange};
for p = 1:2
    subplot(1, 2, p);
    histogram(Parameters(:, p), 500, 'Normalization', 'pdf', 'EdgeColor', 'none');
    xlim(Limits{p});
    title(Titles{p}, 'Interpreter', 'latex', 'FontSize', 13);
    set(gca, 'FontSize', 11, 'FontWeight', 'bold');
end
exportgraphics(Figure, fullfile(FigureFolder, 'paramhist.pdf'), 'ContentType', 'vector');
close(Figure);

% ===== Figure 3: fixed states x_mi and x_hi (CAT bond)
Figure = figure('Visible', 'off', 'Position', [100 100 1200 480]);
Panels = ["mi", "hi"];
PanelTitles = {'period n MSE for center state-x_{mi}', 'period n MSE for extreme state-x_{hi}'};
for p = 1:2
    subplot(1, 2, p);
    plot_panel(CatBond, Panels(p), Style, 5);
    ylim([3e-6 2e-2]);
    title(PanelTitles{p}, 'FontSize', 13);
end
exportgraphics(Figure, fullfile(FigureFolder, 'sig2_lomihi.pdf'), 'ContentType', 'vector');
close(Figure);

% ===== Figure 4: current state (CAT bond)
Figure = figure('Visible', 'off', 'Position', [100 100 900 480]);
plot_panel(CatBond, "current", Style, 5);
title('period n MSE for the current state', 'FontSize', 13);
exportgraphics(Figure, fullfile(FigureFolder, 'sig2n.pdf'), 'ContentType', 'vector');
close(Figure);

% ===== Figure 5: current state (credit risk), with error bars of one
% standard error of the estimated MSE, as in the published figure
Figure = figure('Visible', 'off', 'Position', [100 100 900 480]);
hold on
Lines = gobjects(4, 1);
for e = 1:4
    Rows = Credit(Credit.Estimator == Style.Estimators(e), :);
    Lines(e) = errorbar(Rows.Period, Rows.MSE, Rows.SE, Style.Lines{e}, ...
        'Color', Style.Colours(e, :), 'LineWidth', 3);
end
hold off
set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 13, 'FontWeight', 'bold');
xlim([1 52]); ylim([1e-6 1e-3]);
legend(Lines, Style.Estimators(1:4), 'Location', 'southwest', 'Box', 'off', 'FontSize', 14);
xlabel('Number of periods, n');
title('period n variance for the current state', 'FontSize', 13);
exportgraphics(Figure, fullfile(FigureFolder, 'sig2n_example2.pdf'), 'ContentType', 'vector');
close(Figure);
end


function plot_panel(Summary, Target, Style, NumEstimators)
% Log-log plot of the MSEs of the estimators at one target, with the line
% of slope -1 through the SMC MSE at n = 1 (SMC with n r replications)
hold on
Lines = gobjects(NumEstimators + 1, 1);
for e = 1:NumEstimators
    Rows = Summary(Summary.Target == Target & Summary.Estimator == Style.Estimators(e), :);
    Lines(e) = plot(Rows.Period, Rows.MSE, Style.Lines{e}, 'Color', Style.Colours(e, :), 'LineWidth', 4);
end
SMCFirst = Summary.MSE(Summary.Target == Target & Summary.Estimator == "SMC" & Summary.Period == 1);
Reference = logspace(0, 2, 20);
Lines(end) = plot(Reference, SMCFirst ./ Reference, '-', 'Color', 'k', 'LineWidth', 4);
hold off
set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 13, 'FontWeight', 'bold');
xlim([1 100]);
legend(Lines, [Style.Estimators(1:NumEstimators), "Slope -1"], 'Location', 'southwest', ...
    'Box', 'off', 'FontSize', 14);
xlabel('Number of periods, n');
end
