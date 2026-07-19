/* TEST 00: confirm that the macro library was compiled */

%put NOTE: BLRM_MAIN_EXISTS=%sysmacexist(blrm_main);
%put NOTE: BLRM_MAIN_MC_EXISTS=%sysmacexist(blrm_main_mc);
%put NOTE: BLRM_STAT_EXISTS=%sysmacexist(blrm_stat);
%put NOTE: BLRM_EWOC_EXISTS=%sysmacexist(blrm_ewoc);
%put NOTE: BLRM_SIM_ONE_EXISTS=%sysmacexist(blrm_sim_one);
%put NOTE: BLRM_SIM_N_EXISTS=%sysmacexist(blrm_sim_n);
%put NOTE: BLRM_SIM_SUM_EXISTS=%sysmacexist(blrm_sim_sum);

/*======================================================================*/
/* Test 01: Basic single-chain execution of blrm_main                    */
/*                                                                      */
/* Prerequisite: run sas/blrm_m2.sas in the current SAS session first.  */
/*======================================================================*/

options mprint mlogic;
ods graphics off;


/* Remove outputs from any previous run. */
proc datasets library=work nolist;
    delete test01_input test01_posterior test01_postsum;
quit;


/* Illustrative dose-level data. */
data work.test01_input;
    input dose n dltn;
    datalines;
1 3 0
5 6 1
;
run;


/* Basic single-chain BLRM analysis. */
%blrm_main(
    datain=work.test01_input,
    dataout=work.test01_posterior,
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
    diagnostics=0,
    diagplots=0,
    init_random=0,
    postsumout=work.test01_postsum
);


/* Automated validation. */
%macro check_test01;

    %local nobs miss_alpha miss_beta;

    %if not %sysfunc(exist(work.test01_posterior)) %then %do;
        %put ERROR: TEST01_FAIL - posterior dataset was not created.;
        %return;
    %end;

    proc sql noprint;
        select count(*),
               sum(missing(log_alpha)),
               sum(missing(log_beta))
        into :nobs trimmed,
             :miss_alpha trimmed,
             :miss_beta trimmed
        from work.test01_posterior;
    quit;

    %if &nobs. = 2000
        and &miss_alpha. = 0
        and &miss_beta. = 0
    %then %do;
        %put NOTE: TEST01_PASS - posterior draws are complete.;
    %end;
    %else %do;
        %put ERROR: TEST01_FAIL - nobs=&nobs.
            miss_alpha=&miss_alpha. miss_beta=&miss_beta.;
    %end;

    %if %sysfunc(exist(work.test01_postsum)) %then
        %put NOTE: TEST01_PASS - posterior summary was created.;
    %else
        %put ERROR: TEST01_FAIL - posterior summary was not created.;

%mend check_test01;

%check_test01;


/* Display key outputs. */
proc print data=work.test01_postsum;
    title "Test 01: Posterior Summary";
run;

proc print data=work.test01_posterior(obs=10);
    var Iteration log_alpha log_beta;
    title "Test 01: First 10 Posterior Draws";
run;

title;