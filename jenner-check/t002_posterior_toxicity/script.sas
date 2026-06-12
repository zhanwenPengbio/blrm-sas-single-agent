/* Posterior dose-toxicity summary using the single-agent BLRM transform from
   sas/blrm_m2.sas. The expression

       logistic(log_alpha + exp(log_beta) * log(dose / ref_dose))

   is the model's posterior toxicity probability, written exactly as it
   appears in the repository: %blrm_main (delta / pi_d, lines 42-44),
   %blrm_stat (line 87) and %blrm_ewoc (line 159). Here it is applied to the
   posterior draws in out_mcmc for four candidate doses, then summarized with
   PROC MEANS to give the posterior mean, standard deviation and quartiles of
   the DLT probability at each dose -- the posterior-toxicity summary the
   library is built to produce. */

%let ref_dose = 20;

data post_tox;
  set out_mcmc;
  pi_5  = logistic(log_alpha + exp(log_beta) * log(5  / &ref_dose.));
  pi_10 = logistic(log_alpha + exp(log_beta) * log(10 / &ref_dose.));
  pi_20 = logistic(log_alpha + exp(log_beta) * log(20 / &ref_dose.));
  pi_40 = logistic(log_alpha + exp(log_beta) * log(40 / &ref_dose.));
run;

proc means data=post_tox n mean std p25 median p75 maxdec=4;
  var pi_5 pi_10 pi_20 pi_40;
run;
