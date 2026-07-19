/*======================================================================*/
/* Test 03: Multi-chain BLRM analysis and R-hat                         */
/*                                                                      */
/* Prerequisite: run sas/blrm_m2.sas in the current SAS session first.  */
/*======================================================================*/

options mprint mlogic;
ods graphics off;


/* Remove outputs from any previous run. */
proc datasets library=work nolist nowarn;
    delete test03_input
           test03_posterior
           test03_postsum
           test03_rhat
           test03_chain_diag_ess
           test03_chain_diag_mcse
           test03_chain_diag_autocorr
           test03_chain_diag_geweke
           test03_chain_counts;
quit;


/* Illustrative dose-level data. */
data work.test03_input;
    input dose n dltn;
    datalines;
1 3 0
5 6 1
;
run;


/* Four independent chains. */
%blrm_main_mc(
    datain=work.test03_input,
    dataout=work.test03_posterior,
    ref_dose=5,
    mu1=-1.3863,
    mu2=0,
    v1=4,
    v2=1,
    rho=0,
    alpha=0.05,
    nbi=500,
    nmc=2000,
    nchain=4,
    simming=0,
    seed=9527,
    diagnostics=1,
    chain_diagout=work.test03_chain_diag,
    diagplots=0,
    ac_lags=%str(1 5 10 50),
    init_random=1,
    rhat_cutoff=1.01,
    postsumout=work.test03_postsum,
    diagout=work.test03_rhat,
    debug=0
);


/* Check whether a data set exists and is nonempty. */
%macro assert_nonempty03(ds=);

    %local dsid nobs rc;

    %if not %sysfunc(exist(&ds.)) %then %do;
        %put ERROR: TEST03_FAIL - &ds. was not created.;
        %return;
    %end;

    %let dsid=%sysfunc(open(&ds.));
    %let nobs=%sysfunc(attrn(&dsid.,NOBS));
    %let rc=%sysfunc(close(&dsid.));

    %if &nobs. > 0 %then
        %put NOTE: TEST03_PASS - &ds. contains &nobs. observations.;
    %else
        %put ERROR: TEST03_FAIL - &ds. is empty.;

%mend assert_nonempty03;


/* Verify expected output data sets. */
%assert_nonempty03(ds=work.test03_posterior);
%assert_nonempty03(ds=work.test03_postsum);
%assert_nonempty03(ds=work.test03_rhat);
%assert_nonempty03(ds=work.test03_chain_diag_ess);
%assert_nonempty03(ds=work.test03_chain_diag_mcse);
%assert_nonempty03(ds=work.test03_chain_diag_autocorr);
%assert_nonempty03(ds=work.test03_chain_diag_geweke);


/* Verify pooled posterior structure. */
/*======================================================================*/
/* Corrected Test 03 validation                                         */
/*======================================================================*/


/* Collect pooled posterior information. */
proc sql noprint;
    select count(*),
           count(distinct chain),
           min(chain_iteration),
           max(chain_iteration),
           count(distinct chain_seed),
           min(chain_seed),
           max(chain_seed)
    into :total_draws trimmed,
         :n_chains trimmed,
         :min_iteration trimmed,
         :max_iteration trimmed,
         :n_seeds trimmed,
         :min_seed trimmed,
         :max_seed trimmed
    from work.test03_posterior;
quit;


/* Count retained draws in each chain. */
proc sql;
    create table work.test03_chain_counts as
    select chain,
           chain_seed,
           count(*) as N_Draws
    from work.test03_posterior
    group by chain, chain_seed
    order by chain;
quit;

proc sql noprint;
    select sum(N_Draws ne 2000)
    into :bad_chain_count trimmed
    from work.test03_chain_counts;
quit;


/* Collect R-hat output information. */
proc sql noprint;
    select count(*),
           sum(missing(Rhat)),
           min(N_Chain),
           max(N_Chain)
    into :rhat_rows trimmed,
         :missing_rhat trimmed,
         :min_rhat_chains trimmed,
         :max_rhat_chains trimmed
    from work.test03_rhat;
quit;


/* All macro conditions must be evaluated inside a macro. */
%macro check_test03;

    %if &total_draws. = 8000
        and &n_chains. = 4
        and &min_iteration. = 1
        and &max_iteration. = 2000
        and &n_seeds. = 4
        and &min_seed. = 9527
        and &max_seed. = 9530
    %then %do;
        %put NOTE: TEST03_PASS - pooled posterior structure is correct.;
    %end;
    %else %do;
        %put ERROR: TEST03_FAIL - pooled posterior structure is incorrect.;
        %put ERROR: total_draws=&total_draws. n_chains=&n_chains.;
        %put ERROR: iterations=&min_iteration.-&max_iteration.;
        %put ERROR: seeds=&min_seed.-&max_seed. n_seeds=&n_seeds.;
    %end;


    %if &bad_chain_count. = 0 %then %do;
        %put NOTE: TEST03_PASS - every chain contains 2000 draws.;
    %end;
    %else %do;
        %put ERROR: TEST03_FAIL - one or more chains have an incorrect size.;
    %end;


    %if &rhat_rows. = 2
        and &missing_rhat. = 0
        and &min_rhat_chains. = 4
        and &max_rhat_chains. = 4
    %then %do;
        %put NOTE: TEST03_PASS - R-hat output structure is correct.;
    %end;
    %else %do;
        %put ERROR: TEST03_FAIL - R-hat output structure is incorrect.;
        %put ERROR: rhat_rows=&rhat_rows. missing_rhat=&missing_rhat.;
        %put ERROR: N_Chain range=&min_rhat_chains.-&max_rhat_chains.;
    %end;

%mend check_test03;

%check_test03;


/* Display key outputs. */
proc print data=work.test03_chain_counts noobs;
    title "Test 03: Draws by Chain";
run;

proc print data=work.test03_rhat noobs;
    var parameter
        N_Chain
        N_Draw_Per_Chain
        Rhat
        Rhat_Cutoff
        Converged
        Rhat_Status;
    title "Test 03: Multi-chain R-hat";
run;

proc print data=work.test03_postsum noobs;
    title "Test 03: Pooled Posterior Summary";
run;

title;