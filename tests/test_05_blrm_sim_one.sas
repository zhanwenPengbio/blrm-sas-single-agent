/*======================================================================*/
/* Test 05: One simulated BLRM trial                                    */
/*                                                                      */
/* Tests:                                                               */
/*   - numeric dose mapping                                              */
/*   - one complete simulated trial                                     */
/*   - final cohort-size adjustment                                     */
/*   - DLT, EWOC, and trial-level output structures                      */
/*                                                                      */
/* Prerequisite: run sas/blrm_m2.sas in the current SAS session first.  */
/*======================================================================*/

options mprint mlogic;
ods graphics off;


/* Remove outputs from previous runs. */
proc datasets library=work nolist nowarn;
    delete test05_truth_map_input
           test05_truth_zero
           test05_map
           test05_ewoc_counts
           dlt_record
           ewoc_record
           sim_record;
quit;


/*--------------------------------------------------------------*/
/* Part A: verify dose mapping with deliberately unsorted truth. */
/*--------------------------------------------------------------*/

data work.test05_truth_map_input;
    input dose true_tox;
    datalines;
10 0.30
1  0.10
5  0.20
;
run;

%make_dose_truth_map(
    dose_list=%str(1,5,10),
    true_tox=work.test05_truth_map_input,
    out=work.test05_map
);


/*--------------------------------------------------------------*/
/* Part B: deterministic zero-toxicity simulation scenario.     */
/* The truth data are again deliberately unsorted.              */
/*--------------------------------------------------------------*/

data work.test05_truth_zero;
    input dose true_tox;
    datalines;
10 0
1  0
5  0
;
run;

%blrm_sim_one(
    dose_list=%str(1,5,10),
    ref_dose=5,
    mu1=-2.1972246,
    mu2=0,
    v1=0.25,
    v2=0.25,
    rho=0,
    alpha=0.05,
    nbi=200,
    nmc=500,
    ud=0.16,
    od=0.33,
    ewoc=0.25,
    true_tox=work.test05_truth_zero,
    start_dose=1,
    max_sample=10,
    cohortn=3,
    skip=0,
    seed=20260719,
    sim_batch=101,
    debug=0
);


/*--------------------------------------------------------------*/
/* Collect validation information.                              */
/*--------------------------------------------------------------*/

proc sql noprint;

    /* Dose-list order must be retained despite unsorted truth data. */
    select count(*)
    into :map_rows trimmed
    from work.test05_map;

    select count(*)
    into :bad_map trimmed
    from work.test05_map
    where (dose_id=1 and (dose ne 1  or true_tox ne 0.10))
       or (dose_id=2 and (dose ne 5  or true_tox ne 0.20))
       or (dose_id=3 and (dose ne 10 or true_tox ne 0.30))
       or dose_id not in (1,2,3);


    /* One row per treated cohort is expected. */
    select count(*),
           sum(n),
           sum(dltn),
           min(n),
           max(n),
           count(distinct enroll_batch),
           min(enroll_batch),
           max(enroll_batch)
    into :dlt_rows trimmed,
         :total_sample trimmed,
         :total_dlt trimmed,
         :min_cohort trimmed,
         :max_cohort trimmed,
         :n_batches trimmed,
         :min_batch trimmed,
         :max_batch trimmed
    from work.dlt_record;


    /* Every analysis batch should contain all three candidate doses. */
    create table work.test05_ewoc_counts as
    select enroll_batch,
           count(*) as N_Doses
    from work.ewoc_record
    group by enroll_batch
    order by enroll_batch;

    select count(*),
           sum(N_Doses ne 3)
    into :ewoc_batches trimmed,
         :bad_ewoc_batches trimmed
    from work.test05_ewoc_counts;


    /* Exactly one trial-level result is expected. */
    select count(*),
           min(sample_used),
           max(sample_used),
           sum(missing(dose_select)),
           min(sim_batch),
           max(sim_batch)
    into :sim_rows trimmed,
         :min_sample_used trimmed,
         :max_sample_used trimmed,
         :missing_selection trimmed,
         :min_sim_batch trimmed,
         :max_sim_batch trimmed
    from work.sim_record;

quit;


/*--------------------------------------------------------------*/
/* Automated checks.                                            */
/*--------------------------------------------------------------*/

%macro check_test05;

    %if &map_rows. = 3 and &bad_map. = 0 %then %do;
        %put NOTE: TEST05_PASS - dose and truth mapping is correct.;
    %end;
    %else %do;
        %put ERROR: TEST05_FAIL - dose and truth mapping is incorrect.;
    %end;


    /*
    max_sample=10 and cohortn=3 should produce cohort sizes:
        3, 3, 3, 1
    */
    %if &dlt_rows. = 4
        and &total_sample. = 10
        and &total_dlt. = 0
        and &min_cohort. = 1
        and &max_cohort. = 3
        and &n_batches. = 4
        and &min_batch. = 1
        and &max_batch. = 4
    %then %do;
        %put NOTE: TEST05_PASS - enrollment and final cohort adjustment are correct.;
    %end;
    %else %do;
        %put ERROR: TEST05_FAIL - enrollment structure is incorrect.;
        %put ERROR: rows=&dlt_rows. sample=&total_sample. dlt=&total_dlt.;
        %put ERROR: cohort_range=&min_cohort.-&max_cohort.;
        %put ERROR: batch_range=&min_batch.-&max_batch. n_batches=&n_batches.;
    %end;


    %if &ewoc_batches. = 4 and &bad_ewoc_batches. = 0 %then %do;
        %put NOTE: TEST05_PASS - every enrollment batch has three EWOC dose rows.;
    %end;
    %else %do;
        %put ERROR: TEST05_FAIL - EWOC batch structure is incorrect.;
    %end;


    %if &sim_rows. = 1
        and &min_sample_used. = 10
        and &max_sample_used. = 10
        and &missing_selection. = 0
        and &min_sim_batch. = 101
        and &max_sim_batch. = 101
    %then %do;
        %put NOTE: TEST05_PASS - final simulated-trial record is correct.;
    %end;
    %else %do;
        %put ERROR: TEST05_FAIL - final simulated-trial record is incorrect.;
    %end;

%mend check_test05;

%check_test05;


/*--------------------------------------------------------------*/
/* Display outputs.                                             */
/*--------------------------------------------------------------*/

proc print data=work.test05_map noobs;
    title "Test 05: Dose and Truth Mapping";
run;

proc print data=work.dlt_record noobs;
    title "Test 05: Cohort-level DLT Records";
run;

proc print data=work.test05_ewoc_counts noobs;
    title "Test 05: EWOC Rows by Enrollment Batch";
run;

proc print data=work.sim_record noobs;
    title "Test 05: Final Trial Record";
run;

title;