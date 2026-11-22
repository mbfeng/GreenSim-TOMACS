function price = CatBondPriceVar(PoissonMean, ExpMean, Trigger, Recovery, nTerms)

Pj = poisspdf(0:1:nTerms, PoissonMean);
Fj = gamcdf(Trigger, 0:1:nTerms, ExpMean);
pr = Pj*Fj';
price = (1-Recovery)^2*(pr-pr^2);   % true CatBond price for state 1

end