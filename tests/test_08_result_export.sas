/*======================================================================*/
/* Test 08: Current BLRM Excel export workflow                          */
/*                                                                      */
/* Tests:                                                               */
/*   - workbook creation                                                 */
/*   - INFO, DLT, EWOC, SIM, and OC sheets                              */
/*   - correct observation counts                                       */
/*   - absence of obsolete DOSE_SELECTION and SIM_RESULTS sheets         */
/*                                                                      */
/* Prerequisite: run the updated sas/blrm_export.sas first.             */
/*======================================================================*/

options mprint mlogic;


/* Remove test data from previous runs. */
proc datasets library=work nolist nowarn;
    delete test08_info
           test08_dlt
           test08_ewoc
           test08_sim
           test08_oc
           test08_sheet_check;
quit;


/*--------------------------------------------------------------*/
/* Create representative current-workflow outputs.              */
/*--------------------------------------------------------------*/

data work.test08_info;
    length Scenario $40 Description $100;
    Scenario='Illustrative';
    Description='BLRM export validation';
run;


data work.test08_dlt;
    input sim_batch enroll_batch dose n dltn;
    datalines;
1 1 1 3 0
1 2 5 3 1
;
run;


data work.test08_ewoc;
    length EWOC $5;
    input sim_batch enroll_batch dose
          UD_Mean TT_Mean OD_Mean EWOC $;
    datalines;
1 1 1 0.90 0.10 0.00 True
1 1 5 0.20 0.60 0.20 True
;
run;


data work.test08_sim;
    input sim_batch dose_select sample_used;
    datalines;
1 5 6
;
run;


data work.test08_oc;
    length Metric $50 Parameters $200;

    Parameters=
        'Acceptable: [0.05, 0.33) | Target: [0.16, 0.33) | Over: [0.33, 1]';

    input Metric :$50. Value;

    datalines;
N_Simulations 1
Prop_Correct_Dose 1
Prop_Acceptable_Dose 1
Prop_Over_Toxic_Dose 0
Prop_No_Dose_Selected 0
Avg_Sample_Used 6
Avg_Patients_Over_Toxic 0
;
run;


/*--------------------------------------------------------------*/
/* Export to a temporary workbook in the SAS WORK directory.    */
/*--------------------------------------------------------------*/

%let test08_path=%sysfunc(pathname(work));
%let test08_workbook=&test08_path.\test08_result_export.xlsx;

%result_export(
    path=&test08_path.,
    file_name=test08_result_export,
    info=work.test08_info,
    dlt_data=work.test08_dlt,
    ewoc_data=work.test08_ewoc,
    sim_data=work.test08_sim,
    oc_data=work.test08_oc
);


/*--------------------------------------------------------------*/
/* Confirm that the physical workbook exists.                   */
/*--------------------------------------------------------------*/

filename t08file "&test08_workbook.";

%let workbook_exists=%sysfunc(fexist(t08file));

filename t08file clear;


/*--------------------------------------------------------------*/
/* Reopen the workbook and inspect its sheets.                  */
/*--------------------------------------------------------------*/

libname t08xlsx xlsx "&test08_workbook.";

%let xlsx_libref_rc=%sysfunc(libref(t08xlsx));

%let info_exists=%sysfunc(exist(t08xlsx.INFO));
%let dlt_exists=%sysfunc(exist(t08xlsx.DLT));
%let ewoc_exists=%sysfunc(exist(t08xlsx.EWOC));
%let sim_exists=%sysfunc(exist(t08xlsx.SIM));
%let oc_exists=%sysfunc(exist(t08xlsx.OC));

%let old_selection_exists=
    %sysfunc(exist(t08xlsx.DOSE_SELECTION));

%let old_results_exists=
    %sysfunc(exist(t08xlsx.SIM_RESULTS));


/* Count observations in every expected sheet. */
proc sql noprint;

    select count(*) into :info_n trimmed
    from t08xlsx.INFO;

    select count(*) into :dlt_n trimmed
    from t08xlsx.DLT;

    select count(*) into :ewoc_n trimmed
    from t08xlsx.EWOC;

    select count(*) into :sim_n trimmed
    from t08xlsx.SIM;

    select count(*) into :oc_n trimmed
    from t08xlsx.OC;

quit;


/* Create a compact sheet-validation table. */
data work.test08_sheet_check;
    length Sheet $20;
    length Exists N_Obs 8;

    Sheet='INFO';
    Exists=&info_exists.;
    N_Obs=&info_n.;
    output;

    Sheet='DLT';
    Exists=&dlt_exists.;
    N_Obs=&dlt_n.;
    output;

    Sheet='EWOC';
    Exists=&ewoc_exists.;
    N_Obs=&ewoc_n.;
    output;

    Sheet='SIM';
    Exists=&sim_exists.;
    N_Obs=&sim_n.;
    output;

    Sheet='OC';
    Exists=&oc_exists.;
    N_Obs=&oc_n.;
    output;

    Sheet='DOSE_SELECTION';
    Exists=&old_selection_exists.;
    N_Obs=.;
    output;

    Sheet='SIM_RESULTS';
    Exists=&old_results_exists.;
    N_Obs=.;
    output;
run;


libname t08xlsx clear;


/*--------------------------------------------------------------*/
/* Automated checks.                                            */
/*--------------------------------------------------------------*/

%macro check_test08;

    %if &workbook_exists. = 1 %then %do;
        %put NOTE: TEST08_PASS - Excel workbook was created.;
    %end;
    %else %do;
        %put ERROR: TEST08_FAIL - Excel workbook was not created.;
    %end;


    %if &xlsx_libref_rc. = 0 %then %do;
        %put NOTE: TEST08_PASS - workbook can be reopened with the XLSX engine.;
    %end;
    %else %do;
        %put ERROR: TEST08_FAIL - workbook cannot be reopened.;
    %end;


    %if &info_exists. = 1
        and &dlt_exists. = 1
        and &ewoc_exists. = 1
        and &sim_exists. = 1
        and &oc_exists. = 1
    %then %do;
        %put NOTE: TEST08_PASS - all five current workflow sheets exist.;
    %end;
    %else %do;
        %put ERROR: TEST08_FAIL - one or more required sheets are missing.;
    %end;


    %if &info_n. = 1
        and &dlt_n. = 2
        and &ewoc_n. = 2
        and &sim_n. = 1
        and &oc_n. = 7
    %then %do;
        %put NOTE: TEST08_PASS - all sheet observation counts are correct.;
    %end;
    %else %do;
        %put ERROR: TEST08_FAIL - one or more sheet counts are incorrect.;
    %end;


    %if &old_selection_exists. = 0
        and &old_results_exists. = 0
    %then %do;
        %put NOTE: TEST08_PASS - obsolete sheets are absent.;
    %end;
    %else %do;
        %put ERROR: TEST08_FAIL - obsolete sheets remain in the workbook.;
    %end;

%mend check_test08;

%check_test08;


/* Display validation results. */
proc print data=work.test08_sheet_check noobs;
    title "Test 08: Excel Sheet Validation";
run;

title;