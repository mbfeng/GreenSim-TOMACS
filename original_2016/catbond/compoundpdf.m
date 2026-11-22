function jointpdf = compoundpdf(N, C, lambda, theta)
    jointpdf = zeros(length(N), 1);
    Ind = (N>0);
    
    jointpdf(Ind) = poisspdf(N(Ind), lambda) .* ...
        gampdf(C(Ind), N(Ind), theta);
    jointpdf(~Ind) = poisspdf(N(~Ind), lambda);

end