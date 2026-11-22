function pdf = CopulaPDF(Vt1, Vt2, mu1, mu2, Port)

U = [logncdf(Vt1, mu1 + Port.mu_t(1), Port.sig_t(1)),...
    logncdf(Vt2, mu2 + Port.mu_t(2), Port.sig_t(2))];
U(U<1e-8) = 0;

pdf = copulapdf('t', U, Port.Cor,Port.df).*...
    lognpdf(Vt1, mu1 + Port.mu_t(1), Port.sig_t(1)).*...
    lognpdf(Vt2, mu2 + Port.mu_t(2), Port.sig_t(2));

end