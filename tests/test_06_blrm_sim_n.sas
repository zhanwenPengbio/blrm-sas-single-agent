/*======================================================================*/
/* Test 06: Repeated simulations with blrm_sim_n                        */
/*                                                                      */
/* Tests:                                                               */
/*   - aggregation of three simulated trials                            */
/*   - DLT and EWOC records by sim_batch                                */
/*   - trial-level simulation records                                   */
/*   - total and per-simulation timing output                           */
/*                                                                      */
/* Prerequisite: run sas/blrm_m2.sas in the current SAS session first.  */
/*======================================================================*/

options mprint mlogic;
ods graphics off;


/* Remove outputs from previous runs. */
proc datasets library=work nolist nowarn;
    delete test06_truth_zero
           test06_dlt_by_sim
           test06_ewoc_by_sim
           test06_time_summary
           dlt_all
           ewoc_all
           sim_all
           time_record
           dlt_record
           ewoc_record
           sim_record;
quit;


/* Deterministic zero-toxicity scenario. */
data work.test06_truth_zero;
    input dose true_tox;
    datalines;
10 0
1  0
5  0
;
run;


/* Run three complete simulated trials. */
%blrm_sim_n(
    dose_list=%str(1,5,10),
    ref_dose=5,
    true_tox=work.test06_truth_zero,
    start_dose=1,
    max_sample=10,
    cohortn=3,
    skip=0,
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
    seed=20260719,
    sim_time=3,
    debug=0,
    time_out=work.test06_time_summary,
    time_sum=0
);


/*--------------------------------------------------------------*/
/* Build per-simulation summaries for validation.               */
/*--------------------------------------------------------------*/

proc sql;

    create table work.test06_dlt_by_sim as
    select sim_batch,
           count(*) as N_Cohorts,
           sum(n) as Total_Sample,
           sum(dltn) as Total_DLT,
           min(n) as Min_Cohort,
           max(n) as Max_Cohort
    from work.dlt_all
    group by sim_batch
    order by sim_batch;


    create table work.test06_ewoc_by_sim as
    select sim_batch,
           count(*) as N_EWOC_Rows,
           count(distinct enroll_batch) as N_Enrollment_Batches
    from work.ewoc_all
    group by sim_batch
    order by sim_batch;

quit;


/*--------------------------------------------------------------*/
/* Collect validation information.                              */
/*--------------------------------------------------------------*/

proc sql noprint;

    select count(*),
           sum(N_Cohorts ne 4),
           sum(Total_Sample ne 10),
           sum(Total_DLT ne 0),
           sum(Min_Cohort ne 1),
           sum(Max_Cohort ne 3)
    into :dlt_sim_rows trimmed,
         :bad_cohort_rows trimmed,
         :bad_sample_totals trimmed,
         :bad_dlt_totals trimmed,
         :bad_min_cohort trimmed,
         :bad_max_cohort trimmed
    from work.test06_dlt_by_sim;


    select count(*),
           sum(N_EWOC_Rows ne 12),
           sum(N_Enrollment_Batches ne 4)
    into :ewoc_sim_rows trimmed,
         :bad_ewoc_rows trimmed,
         :bad_enrollment_batches trimmed
    from work.test06_ewoc_by_sim;


    select count(*),
           count(distinct sim_batch),
           min(sim_batch),
           max(sim_batch),
           sum(sample_used ne 10),
           sum(missing(dose_select))
    into :sim_all_rows trimmed,
         :sim_batch_count trimmed,
         :min_sim_batch trimmed,
         :max_sim_batch trimmed,
         :bad_sample_used trimmed,
         :missing_dose_select trimmed
    from work.sim_all;


    /*
    test06_time_summary should contain:
      - one total-time row with missing sim_batch
      - three per-simulation rows with nonmissing sim_batch
    */
    select count(*),
           sum(missing(sim_batch)),
           sum(not missing(sim_batch)),
           sum(not missing(sim_batch) and missing(elapsed_sec)),
           sum(
               missing(sim_batch)
               and total_elapsed_sec > 0
               and avg_sec_per_sim > 0
           )
    into :time_rows trimmed,
         :total_time_rows trimmed,
         :individual_time_rows trimmed,
         :missing_individual_time trimmed,
         :valid_total_time_rows trimmed
    from work.test06_time_summary;

quit;


/*--------------------------------------------------------------*/
/* Automated checks.                                            */
/*--------------------------------------------------------------*/

%macro check_test06;

    %if &dlt_sim_rows. = 3
        and &bad_cohort_rows. = 0
        and &bad_sample_totals. = 0
        and &bad_dlt_totals. = 0
        and &bad_min_cohort. = 0
        and &bad_max_cohort. = 0
    %then %do;
        %put NOTE: TEST06_PASS - DLT records are correctly aggregated by simulation.;
    %end;
    %else %do;
        %put ERROR: TEST06_FAIL - DLT aggregation is incorrect.;
    %end;


    %if &ewoc_sim_rows. = 3
        and &bad_ewoc_rows. = 0
        and &bad_enrollment_batches. = 0
    %then %do;
        %put NOTE: TEST06_PASS - EWOC records are correctly aggregated by simulation.;
    %end;
    %else %do;
        %put ERROR: TEST06_FAIL - EWOC aggregation is incorrect.;
    %end;


    %if &sim_all_rows. = 3
        and &sim_batch_count. = 3
        and &min_sim_batch. = 1
        and &max_sim_batch. = 3
        and &bad_sample_used. = 0
        and &missing_dose_select. = 0
    %then %do;
        %put NOTE: TEST06_PASS - trial-level simulation records are correct.;
    %end;
    %else %do;
        %put ERROR: TEST06_FAIL - trial-level simulation records are incorrect.;
    %end;


    %if &time_rows. = 4
        and &total_time_rows. = 1
        and &individual_time_rows. = 3
        and &missing_individual_time. = 0
        and &valid_total_time_rows. = 1
    %then %do;
        %put NOTE: TEST06_PASS - timing output is correct.;
    %end;
    %else %do;
        %put ERROR: TEST06_FAIL - timing output is incorrect.;
        %put ERROR: time_rows=&time_rows.
            total_rows=&total_time_rows.
            individual_rows=&individual_time_rows.;
    %end;

%mend check_test06;

%check_test06;


/*--------------------------------------------------------------*/
/* Display outputs.                                             */
/*--------------------------------------------------------------*/

proc print data=work.test06_dlt_by_sim noobs;
    title "Test 06: DLT Records by Simulation";
run;

proc print data=work.test06_ewoc_by_sim noobs;
    title "Test 06: EWOC Records by Simulation";
run;

proc print data=work.sim_all noobs;
    title "Test 06: Trial-level Simulation Results";
run;

proc print data=work.test06_time_summary noobs;
    title "Test 06: Timing Summary";
run;

title;