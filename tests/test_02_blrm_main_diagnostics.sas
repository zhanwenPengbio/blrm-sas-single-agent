/*======================================================================*/
/* Test 02: Single-chain MCMC diagnostics                               */
/*                                                                      */
/* Prerequisite: run sas/blrm_m2.sas in the current SAS session first.  */
/*======================================================================*/

options mprint mlogic;
ods graphics off;


/* Remove outputs from any previous run. */
proc datasets library=work nolist;
    delete test02_input
           test02_posterior
           test02_postsum
           test02_diag_ess
           test02_diag_mcse
           test02_diag_autocorr
           test02_diag_geweke;
quit;


/* Illustrative dose-level data. */
data work.test02_input;
    input dose n dltn;
    datalines;
1 3 0
5 6 1
;
run;


/* Run one chain with diagnostics enabled. */
%blrm_main(
    datain=work.test02_input,
    dataout=work.test02_posterior,
    ref_dose=5,
    mu1=-1.3863,
    mu2=0,
    v1=4,
    v2=1,
    rho=0,
    alpha=0.05,
    nbi=500,
    nmc=2000,
    simming=0,
    seed=9527,
    diagnostics=1,
    diagout=work.test02_diag,
    diagplots=0,
    ac_lags=%str(1 5 10 50),
    init_random=0,
    postsumout=work.test02_postsum
);


/* Check whether a data set exists and contains observations. */
%macro assert_nonempty(ds=);

    %local dsid nobs rc;

    %if not %sysfunc(exist(&ds.)) %then %do;
        %put ERROR: TEST02_FAIL - &ds. was not created.;
        %return;
    %end;

    %let dsid=%sysfunc(open(&ds.));
    %let nobs=%sysfunc(attrn(&dsid.,NOBS));
    %let rc=%sysfunc(close(&dsid.));

    %if &nobs. > 0 %then
        %put NOTE: TEST02_PASS - &ds. contains &nobs. observations.;
    %else
        %put ERROR: TEST02_FAIL - &ds. is empty.;

%mend assert_nonempty;


/* Validate all expected outputs. */
%assert_nonempty(ds=work.test02_posterior);
%assert_nonempty(ds=work.test02_postsum);
%assert_nonempty(ds=work.test02_diag_ess);
%assert_nonempty(ds=work.test02_diag_mcse);
%assert_nonempty(ds=work.test02_diag_autocorr);
%assert_nonempty(ds=work.test02_diag_geweke);


/* Display diagnostic outputs. */
proc print data=work.test02_diag_ess;
    title "Test 02: Effective Sample Size";
run;

proc print data=work.test02_diag_mcse;
    title "Test 02: Monte Carlo Standard Error";
run;

proc print data=work.test02_diag_autocorr(obs=20);
    title "Test 02: Autocorrelation";
run;

proc print data=work.test02_diag_geweke;
    title "Test 02: Geweke Diagnostic";
run;

title;