function test_sk()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test_sk.m                                                                   %
% Checks of the stochastic kriging implementation (sk_fit.m, sk_predict.m):   %
% (1) with tiny intrinsic variances it interpolates the responses;            %
% (2) a surface that is exactly linear is reproduced by the trend;            %
% (3) a warm start from the optimum returns the same fit.                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CodeFolder = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(CodeFolder, 'functions'));
Stream = RandStream('mrg32k3a', 'Seed', 41);
Design = [2 + 2 * rand(Stream, 30, 1), 4 + 2 * rand(Stream, 30, 1)];   % (lambda, theta)
Surface = @(X) 0.9 - 0.05 * (X(:, 1) - 3).^2 + 0.03 * sin(2 * X(:, 2));

% (1) interpolation
Responses = Surface(Design);
Fit = sk_fit(Design, Responses, 1e-12 * ones(30, 1));
Interpolation = max(abs(sk_predict(Fit, Design) - Responses));
fprintf('test_sk: interpolation error %.2e\n', Interpolation);
assert(Interpolation < 1e-4, 'sk_fit does not interpolate noiseless responses');

% (2) linear surface
Linear = @(X) 0.2 + 0.1 * X(:, 1) - 0.05 * X(:, 2);
Fit = sk_fit(Design, Linear(Design), 1e-6 * ones(30, 1));
NewPoints = [2.5 4.5; 3.7 5.9; 3 5];
LinearError = max(abs(sk_predict(Fit, NewPoints) - Linear(NewPoints)));
fprintf('test_sk: linear surface error %.2e\n', LinearError);
assert(LinearError < 1e-4, 'sk_fit does not reproduce a linear surface');

% (3) warm start at the optimum
Noisy = Surface(Design) + 0.01 * randn(Stream, 30, 1);
Variances = 1e-4 * ones(30, 1);
Cold = sk_fit(Design, Noisy, Variances);
Warm = sk_fit(Design, Noisy, Variances, Cold);
WarmDifference = max(abs(sk_predict(Cold, NewPoints) - sk_predict(Warm, NewPoints)));
fprintf('test_sk: warm start difference %.2e\n', WarmDifference);
assert(WarmDifference < 1e-6, 'a warm start at the optimum changes the fit');
end
