options obs=100;

/* Mock posterior draws of the BLRM log-parameters, shaped like the OUTPOST=
   data set that %blrm_main (PROC MCMC) produces in sas/blrm_m2.sas: two
   columns, log_alpha and log_beta. Drawn here from a small Normal cloud so
   the bundle is self-contained and deterministic (fixed seed) rather than
   requiring an MCMC fit. */
data out_mcmc;
  call streaminit(20240611);
  do _i = 1 to 60;
    log_alpha = -0.85 + 0.30 * rand('normal');
    log_beta  =  0.00 + 0.20 * rand('normal');
    output;
  end;
  keep log_alpha log_beta;
run;
