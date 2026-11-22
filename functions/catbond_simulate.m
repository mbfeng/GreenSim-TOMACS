function [Losses, Payoffs] = catbond_simulate(Model, Parameters, NumReplications, Stream)
% Simulation of the CAT bond (Algorithm 3): r replications of the insured
% loss and of the bond's payoff per dollar, given the parameters of a state
% Input(s):
% -- Model: struct from catbond_settings.m
% -- Parameters: 1 x 2, (lambda, theta)
% -- NumReplications: r
% -- Stream: RandStream
% Output(s):
% -- Losses: r x 1, insured losses Y = sum_{i=1}^M Z_i, M ~ Poisson(lambda),
%    Z_i ~ exponential with mean theta; Y = 0 if M = 0
% -- Payoffs: r x 1, F(Y) = 1{Y <= K} + p 1{Y > K}
% The number of claims M and, given M, the loss Y ~ Gamma(M, theta) are drawn
% by inversion from the given stream (poissrnd and gamrnd use the global
% stream), so each macro replication has its own reproducible stream.

Lambda = Parameters(1); Theta = Parameters(2);
Counts = poissinv(rand(Stream, NumReplications, 1), Lambda);   % M
Uniforms = rand(Stream, NumReplications, 1);                   % for the losses
Losses = zeros(NumReplications, 1);
Positive = Counts > 0;
Losses(Positive) = gaminv(Uniforms(Positive), Counts(Positive), Theta); % Y | M ~ Gamma(M, theta)
Payoffs = catbond_payoff(Model, Losses);
end
