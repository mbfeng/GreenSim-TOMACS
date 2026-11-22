function Parameters = catbond_parameters(Model, States)
% Sigmoidal transformation psi from states to loss-model parameters
% (Section 6.1)
% Input(s):
% -- Model: struct from catbond_settings.m
% -- States: N x 2 states X_n
% Output(s):
% -- Parameters: N x 2, (lambda, theta) = psi(X)
% formula: psi(x) = lower + (upper - lower) / (1 + exp(-x)), per coordinate

Lower = [Model.LambdaRange(1), Model.ThetaRange(1)];
Upper = [Model.LambdaRange(2), Model.ThetaRange(2)];
Parameters = Lower + (Upper - Lower) ./ (1 + exp(-States));
end
