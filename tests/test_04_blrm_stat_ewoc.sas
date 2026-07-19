/*======================================================================*/
/* Test 04: blrm_stat and blrm_ewoc                                     */
/*                                                                      */
/* Uses deterministic posterior draws with known expected results.      */
/* Prerequisite: run sas/blrm_m2.sas in the current SAS session first.  */
/*======================================================================*/

options mprint mlogic;


/* Remove outputs from any previous run. */
proc datasets library=work nolist nowarn;
    delete test04_draws
           test04_stat
           test04_ewoc;
quit;


/*
At the reference dose, the four posterior toxicity draws are:

    0.10, 0.20, 0.30, 0.40

log_beta=0 implies exp(log_beta)=1.
*/
data work.test04_draws;
    do pi_ref = 0.10, 0.20, 0.30, 0.40;
        Iteration + 1;
        log_alpha = log(pi_ref / (1 - pi_ref));
        log_beta = 0;
        output;
    end;

    drop pi_ref;
run;


/* Posterior toxicity summaries at three doses. */
%blrm_stat(
    dose_list=%str(1,5,10),
    ref_dose=5,
    mcmcout=work.test04_draws,
    dataout=work.test04_stat,
    debug=0
);


/* EWOC summaries at the same doses. */
%blrm_ewoc(
    dose_list=%str(1,5,10),
    ref_dose=5,
    mcmcout=work.test04_draws,
    dataout=work.test04_ewoc,
    ud=0.16,
    od=0.33,
    ewoc=0.25,
    debug=0
);


/* Collect results for automated validation. */
proc sql noprint;

    select count(*),
           sum(missing(Mean))
    into :stat_rows trimmed,
         :stat_missing trimmed
    from work.test04_stat;

    select Mean
    into :mean_dose5 trimmed
    from work.test04_stat
    where dose = 5;

    select count(*),
           sum(abs(UD_Mean + TT_Mean + OD_Mean - 1) > 1e-12)
    into :ewoc_rows trimmed,
         :bad_probability_sum trimmed
    from work.test04_ewoc;

    select UD_Mean,
           TT_Mean,
           OD_Mean,
           strip(EWOC)
    into :ud1 trimmed,
         :tt1 trimmed,
         :od1 trimmed,
         :flag1 trimmed
    from work.test04_ewoc
    where dose = 1;

    select UD_Mean,
           TT_Mean,
           OD_Mean,
           strip(EWOC)
    into :ud5 trimmed,
         :tt5 trimmed,
         :od5 trimmed,
         :flag5 trimmed
    from work.test04_ewoc
    where dose = 5;

    select UD_Mean,
           TT_Mean,
           OD_Mean,
           strip(EWOC)
    into :ud10 trimmed,
         :tt10 trimmed,
         :od10 trimmed,
         :flag10 trimmed
    from work.test04_ewoc
    where dose = 10;

quit;


/* Automated checks. */
%macro check_test04;

    %if &stat_rows. = 3 and &stat_missing. = 0 %then %do;
        %put NOTE: TEST04_PASS - blrm_stat returned three complete dose rows.;
    %end;
    %else %do;
        %put ERROR: TEST04_FAIL - blrm_stat output is incomplete.;
    %end;


    %if %sysevalf(
        &mean_dose5. > 0.249999 and
        &mean_dose5. < 0.250001
    ) %then %do;
        %put NOTE: TEST04_PASS - posterior mean at reference dose is 0.25.;
    %end;
    %else %do;
        %put ERROR: TEST04_FAIL - reference-dose mean=&mean_dose5.;
    %end;


    %if &ewoc_rows. = 3 and &bad_probability_sum. = 0 %then %do;
        %put NOTE: TEST04_PASS - EWOC probabilities sum to one at every dose.;
    %end;
    %else %do;
        %put ERROR: TEST04_FAIL - EWOC probability output is invalid.;
    %end;


    /*
    Expected result at dose 1:
        UD=1, TT=0, OD=0, EWOC=True
    */
    %if &ud1. = 1
        and &tt1. = 0
        and &od1. = 0
        and %upcase(&flag1.) = TRUE
    %then %do;
        %put NOTE: TEST04_PASS - dose 1 EWOC classification is correct.;
    %end;
    %else %do;
        %put ERROR: TEST04_FAIL - dose 1 EWOC classification is incorrect.;
    %end;


    /*
    Expected result at reference dose 5:
        UD=0.25, TT=0.50, OD=0.25, EWOC=True
    */
    %if &ud5. = 0.25
        and &tt5. = 0.5
        and &od5. = 0.25
        and %upcase(&flag5.) = TRUE
    %then %do;
        %put NOTE: TEST04_PASS - dose 5 EWOC classification is correct.;
    %end;
    %else %do;
        %put ERROR: TEST04_FAIL - dose 5 EWOC classification is incorrect.;
    %end;


    /*
    Expected result at dose 10:
        UD=0, TT=0.25, OD=0.75, EWOC=False
    */
    %if &ud10. = 0
        and &tt10. = 0.25
        and &od10. = 0.75
        and %upcase(&flag10.) = FALSE
    %then %do;
        %put NOTE: TEST04_PASS - dose 10 EWOC classification is correct.;
    %end;
    %else %do;
        %put ERROR: TEST04_FAIL - dose 10 EWOC classification is incorrect.;
    %end;

%mend check_test04;

%check_test04;


/* Display outputs. */
proc print data=work.test04_stat noobs;
    title "Test 04: Posterior Toxicity Summary";
run;

proc print data=work.test04_ewoc noobs;
    title "Test 04: EWOC Summary";
run;

title;