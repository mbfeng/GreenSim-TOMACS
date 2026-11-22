function Payoffs = catbond_payoff(Model, Losses)
% Payoff of the CAT bond per dollar of face value (Section 6.1)
% Input(s):
% -- Model: struct from catbond_settings.m
% -- Losses: insured losses Y
% Output(s):
% -- Payoffs: F(Y) = 1{Y <= K} + p 1{Y > K}
Payoffs = Model.Recovery + (1 - Model.Recovery) * (Losses <= Model.Trigger);
end
