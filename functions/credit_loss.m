function Loss = credit_loss(Model, Scenarios)
% Portfolio loss L_1 + L_2 at the risk horizon t (Section 6.2)
% Input(s):
% -- Model: struct from credit_settings.m
% -- Scenarios: N x 2 asset values S_t
% Output(s):
% -- Loss: N x 1, sum_i ell_i(S_{t,i}),
%    ell_i(s) = a_i if s < D_i, else a_i e^{-r(T-t)} Phi(-d_i),
%    d_i = (ln(s/D_i) + (r - varsigma_i^2/2)(T - t)) / (varsigma_i sqrt(T - t))
% (the risk-neutral price of the loss at maturity T if company i defaults
% then; Phi(-d_i) is the risk-neutral probability of default at T)

Horizon = Model.Evaluation;                                   % T - t
Default = Scenarios < Model.Debt;                             % 1{S_t,i < D_i}
D2 = (log(Scenarios ./ Model.Debt) + (Model.RiskFree - Model.Volatility.^2 / 2) * Horizon) ...
    ./ (Model.Volatility * sqrt(Horizon));
DigitalPut = exp(-Model.RiskFree * Horizon) * normcdf(-D2);
Loss = sum(max(Default, DigitalPut) .* Model.LossGivenDefault, 2);
end
