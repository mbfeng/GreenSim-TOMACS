function jointpdf = besselpdf(C, lambda, theta)
% marginal pdf of Poisson-Exponetial compound loss using 
% modified bessel function of the first kind

    jointpdf = zeros(length(C), 1);
    Ind = (C>0);
    
    TempArg = sqrt(lambda*C(Ind)/theta);
    
    jointpdf(Ind) = (lambda*besseli(1,2*TempArg))./...
        (theta*exp(lambda+C(Ind)/theta).*TempArg);
    jointpdf(~Ind) = poisspdf(0, lambda);

end