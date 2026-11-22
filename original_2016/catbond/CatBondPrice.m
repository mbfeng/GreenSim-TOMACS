function price = CatBondPrice(PoissonMean, ExpMean, Trigger, Recovery, nTerms)

Pj = poisspdf(0:1:nTerms, PoissonMean);
Fj = gamcdf(Trigger, 0:1:nTerms, ExpMean);
price = Recovery + (1-Recovery)*(Pj*Fj');   % true CatBond price for state 1

end