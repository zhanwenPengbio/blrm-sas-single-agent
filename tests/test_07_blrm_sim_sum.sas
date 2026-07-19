/*======================================================================*/
/* Test 07: Operating-characteristic summary with blrm_sim_sum          */
/*                                                                      */
/* Tests:                                                               */
/*   - all simulated trials remain in the denominator                    */
/*   - missing dose selection is counted as zero, not missing            */
/*   - target, acceptable, and over-toxic selection probabilities        */
/*   - trials without DLT records receive zero over-toxic exposure       */
/*                                                                      */
/* Prerequisite: run sas/blrm_m2.sas in the current SAS session first.  */
/*======================================================================*/

options mprint mlogic;


/* Remove outputs from previous runs. */
proc datasets library=work nolist nowarn;
    delete test07_truth
           test07_sim_data
           test07_dlt_data
           test07_oc;
quit;


/*--------------------------------------------------------------*/
/* True toxicity categories                                     */
/*                                                              */
/* Dose 1:  below acceptable range                              */
/* Dose 2:  acceptable, but below target range                  */
/* Dose 5:  target and acceptable                               */
/* Dose 10: over-toxic                                          */
/*--------------------------------------------------------------*/

data work.test07_truth;
    input dose true_tox;
    datalines;
10 0.40
1  0.02
5  0.20
2  0.10
;
run;


/*--------------------------------------------------------------*/
/* Four simulated trials                                        */
/*                                                              */
/* Trial 1 selects target dose 5                                */
/* Trial 2 selects acceptable dose 2                            */
/* Trial 3 selects over-toxic dose 10                           */
/* Trial 4 selects no dose                                      */
/*--------------------------------------------------------------*/

data work.test07_sim_data;
    input sim_batch dose_select sample_used;
    datalines;
1 5  10
2 2   9
3 10  8
4 .   6
;
run;


/*--------------------------------------------------------------*/
/* Cohort-level exposure records                                */
/*                                                              */
/* Trial 3 exposes 5 patients to over-toxic dose 10.            */
/* Trial 4 deliberately has no DLT record.                      */
/*--------------------------------------------------------------*/

data work.test07_dlt_data;
    input sim_batch dose n;
    datalines;
1 1  3
1 5  7
2 1  3
2 2  6
3 5  3
3 10 5
;
run;


/* Produce operating-characteristic summary. */
%blrm_sim_sum(
    sim_data=work.test07_sim_data,
    dlt_data=work.test07_dlt_data,
    true_tox=work.test07_truth,
    target_range_low=0.16,
    target_range_high=0.33,
    acceptable_low=0.05,
    out_summary=work.test07_oc,
    debug=0
);


/*--------------------------------------------------------------*/
/* Collect OC results.                                          */
/*--------------------------------------------------------------*/

proc sql noprint;

    select count(*),
           count(distinct Metric),
           sum(missing(Value))
    into :oc_rows trimmed,
         :oc_unique_metrics trimmed,
         :oc_missing_values trimmed
    from work.test07_oc;


    select Value into :n_simulations trimmed
    from work.test07_oc
    where Metric='N_Simulations';


    select Value into :prop_correct trimmed
    from work.test07_oc
    where Metric='Prop_Correct_Dose';


    select Value into :prop_acceptable trimmed
    from work.test07_oc
    where Metric='Prop_Acceptable_Dose';


    select Value into :prop_over_toxic trimmed
    from work.test07_oc
    where Metric='Prop_Over_Toxic_Dose';


    select Value into :prop_no_selection trimmed
    from work.test07_oc
    where Metric='Prop_No_Dose_Selected';


    select Value into :avg_sample_used trimmed
    from work.test07_oc
    where Metric='Avg_Sample_Used';


    select Value into :avg_over_toxic_exposure trimmed
    from work.test07_oc
    where Metric='Avg_Patients_Over_Toxic';

quit;


/*--------------------------------------------------------------*/
/* Automated checks.                                            */
/*--------------------------------------------------------------*/

%macro check_test07;

    %if &oc_rows. = 7
        and &oc_unique_metrics. = 7
        and &oc_missing_values. = 0
    %then %do;
        %put NOTE: TEST07_PASS - all seven OC metrics were created.;
    %end;
    %else %do;
        %put ERROR: TEST07_FAIL - OC metric structure is incorrect.;
    %end;


    /*
    One of four trials selected the target dose:
        1 / 4 = 0.25
    */
    %if &n_simulations. = 4
        and &prop_correct. = 0.25
    %then %do;
        %put NOTE: TEST07_PASS - target-dose selection probability is correct.;
    %end;
    %else %do;
        %put ERROR: TEST07_FAIL - target-dose selection probability is incorrect.;
    %end;


    /*
    Two of four trials selected acceptable doses:
        dose 5 and dose 2
        2 / 4 = 0.50
    */
    %if &prop_acceptable. = 0.5 %then %do;
        %put NOTE: TEST07_PASS - acceptable-dose selection probability is correct.;
    %end;
    %else %do;
        %put ERROR: TEST07_FAIL - acceptable-dose probability is incorrect.;
    %end;


    /*
    One over-toxic selection and one missing selection:
        1 / 4 = 0.25 for each
    */
    %if &prop_over_toxic. = 0.25
        and &prop_no_selection. = 0.25
    %then %do;
        %put NOTE: TEST07_PASS - over-toxic and no-selection probabilities are correct.;
    %end;
    %else %do;
        %put ERROR: TEST07_FAIL - over-toxic or no-selection probability is incorrect.;
    %end;


    /*
    Average sample size:
        (10 + 9 + 8 + 6) / 4 = 8.25
    */
    %if &avg_sample_used. = 8.25 %then %do;
        %put NOTE: TEST07_PASS - average sample size is correct.;
    %end;
    %else %do;
        %put ERROR: TEST07_FAIL - average sample size is incorrect.;
    %end;


    /*
    Only trial 3 has over-toxic exposure:
        (0 + 0 + 5 + 0) / 4 = 1.25

    Trial 4 has no DLT records but must remain in the denominator.
    */
    %if &avg_over_toxic_exposure. = 1.25 %then %do;
        %put NOTE: TEST07_PASS - over-toxic patient exposure is correct.;
    %end;
    %else %do;
        %put ERROR: TEST07_FAIL - over-toxic patient exposure is incorrect.;
    %end;

%mend check_test07;

%check_test07;


/* Display final OC output. */
proc print data=work.test07_oc noobs;
    var Metric Value Parameters;
    title "Test 07: Operating-Characteristic Summary";
run;

title;