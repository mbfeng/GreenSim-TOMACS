% genpaths.m
% Generate sample paths for later simulations

% param;
PathsMu = norminv(copularnd('t', Port.Cor, Port.df, NPeriods * NPaths));
TrueProb = zeros(NPeriods * NPaths, 1);

for ind = 1:(NPeriods * NPaths)
%     tic
    
    if(mod(ind,NPeriods) == 1)
        PathsMu(ind,:) = log(Port.S0);
        ind/NPeriods
    else
        PathsMu(ind,:) = PathsMu(ind-1,:) + Port.mu_delt + Port.sig_delt.*PathsMu(ind,:);
    end
    TrueProb(ind) = TrueProbFunc(exp(PathsMu(ind,1)), exp(PathsMu(ind,2)));
    
%     toc
end

% clearvars -except Port ExParam PathsMu TrueProb Mu_t Sig_t
save('TrueProb_Apr2016.mat')

