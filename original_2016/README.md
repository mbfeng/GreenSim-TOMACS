# Original code (2016)

The MATLAB code of 2016 that produced the published results, kept for reference (not maintained;
the code in the parent folder replaces it).

- `catbond/`: CAT bond example (Section 6.1). `File1_Initialization.m` (parameters; note that the
  saved file sets the trigger `K = 16`, whereas the runs behind the published figures used
  `K = 25`), `File2_SimPaths.m` (paths and true prices), `File3_SimModel_CRN.m` (SMC, ILR, EWLR
  (`WLR`), MLR and GSK (`SK1`) with common random numbers), `File3_SimModel_CRN_EB.m` (a later
  edit that recomputes SMC and an empirical Bayes shrinkage estimator not used in the paper),
  `File4_Plots.m`, and the price, variance and density functions. `File0_CATBond_main.m` is the
  combined script of the code submitted with the paper.
- `credit/`: credit risk example (Section 6.2). `param.m`, `genpaths.m` (paths and true
  probabilities), `MC_MIS_LimitedMem.m` (SMC, ILR, EWLR, MLR), `MC_MIS_LimitedMem_EB.m` (later edit,
  as above), `File4_Plots.m`, `CopulaPDF.m`, `LossAtt.m`; `File0_CreditRisk_main.m` as submitted.

Not included: the stochastic kriging package called by these scripts (`SKfit`, `SKpredict`,
`SK_GIBF_*`; © 2009 B. L. Nelson, J. Staum, E. Baysal and W. Xie, distributed for personal use
only) and `tight_subplot.m` (MATLAB File Exchange). The parent folder has an independent
implementation of stochastic kriging (`functions/sk_fit.m`, `functions/sk_predict.m`).
