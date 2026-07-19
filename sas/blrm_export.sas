/*======================================================================*/
/* result_export                                                        */
/*                                                                      */
/* path      : output directory                                         */
/* file_name : workbook name without the .xlsx extension                */
/* info      : analysis or scenario information data set                */
/* dlt_data  : cohort-level DLT records                                 */
/* ewoc_data : cohort-level EWOC records                                */
/* sim_data  : trial-level simulation records                           */
/* oc_data   : long-format OC summary from blrm_sim_sum                  */
/*                                                                      */
/* The macro validates every input before creating the workbook.        */
/* If a workbook with the same name already exists, it is replaced.     */
/*======================================================================*/

%macro result_export(
    path=,
    file_name=,
    info=,
    dlt_data=dlt_all,
    ewoc_data=ewoc_all,
    sim_data=sim_all,
    oc_data=
);

    %local _outfile _delete_rc;


    /*--------------------------------------------------------------*/
    /* Validate required parameters.                                */
    /*--------------------------------------------------------------*/

    %if %length(%superq(path)) = 0 %then %do;
        %put ERROR: result_export: path must be specified.;
        %return;
    %end;

    %if %length(%superq(file_name)) = 0 %then %do;
        %put ERROR: result_export: file_name must be specified.;
        %return;
    %end;

    %if %length(%superq(info)) = 0 %then %do;
        %put ERROR: result_export: info must be specified.;
        %return;
    %end;

    %if %length(%superq(oc_data)) = 0 %then %do;
        %put ERROR: result_export: oc_data must be specified.;
        %return;
    %end;


    /*--------------------------------------------------------------*/
    /* Validate all input data sets before writing any sheet.       */
    /*--------------------------------------------------------------*/

    %if not %sysfunc(exist(&info.)) %then %do;
        %put ERROR: result_export: info data set &info. does not exist.;
        %return;
    %end;

    %if not %sysfunc(exist(&dlt_data.)) %then %do;
        %put ERROR: result_export: dlt_data &dlt_data. does not exist.;
        %return;
    %end;

    %if not %sysfunc(exist(&ewoc_data.)) %then %do;
        %put ERROR: result_export: ewoc_data &ewoc_data. does not exist.;
        %return;
    %end;

    %if not %sysfunc(exist(&sim_data.)) %then %do;
        %put ERROR: result_export: sim_data &sim_data. does not exist.;
        %return;
    %end;

    %if not %sysfunc(exist(&oc_data.)) %then %do;
        %put ERROR: result_export: oc_data &oc_data. does not exist.;
        %return;
    %end;


    /*--------------------------------------------------------------*/
    /* Construct output filename and remove any existing workbook.  */
    /*--------------------------------------------------------------*/

    %let _outfile=&path.\&file_name..xlsx;
    %let _delete_rc=0;

    filename _blrmxp "&_outfile.";

    data _null_;
        if fexist('_blrmxp') then do;
            rc=fdelete('_blrmxp');
            call symputx('_delete_rc', rc, 'L');
        end;
        else do;
            call symputx('_delete_rc', 0, 'L');
        end;
    run;

    filename _blrmxp clear;

    %if &_delete_rc. ne 0 %then %do;
        %put ERROR: result_export: existing workbook could not be replaced.;
        %put ERROR: result_export: &_outfile.;
        %return;
    %end;


    /*--------------------------------------------------------------*/
    /* Export current BLRM workflow outputs.                         */
    /*--------------------------------------------------------------*/

    proc export data=&info.
        outfile="&_outfile."
        dbms=xlsx
        replace;
        sheet="INFO";
    run;

    %if &syserr. > 4 %then %do;
        %put ERROR: result_export: failed while writing INFO.;
        %return;
    %end;


    proc export data=&dlt_data.
        outfile="&_outfile."
        dbms=xlsx
        replace;
        sheet="DLT";
    run;

    %if &syserr. > 4 %then %do;
        %put ERROR: result_export: failed while writing DLT.;
        %return;
    %end;


    proc export data=&ewoc_data.
        outfile="&_outfile."
        dbms=xlsx
        replace;
        sheet="EWOC";
    run;

    %if &syserr. > 4 %then %do;
        %put ERROR: result_export: failed while writing EWOC.;
        %return;
    %end;


    proc export data=&sim_data.
        outfile="&_outfile."
        dbms=xlsx
        replace;
        sheet="SIM";
    run;

    %if &syserr. > 4 %then %do;
        %put ERROR: result_export: failed while writing SIM.;
        %return;
    %end;


    proc export data=&oc_data.
        outfile="&_outfile."
        dbms=xlsx
        replace;
        sheet="OC";
    run;

    %if &syserr. > 4 %then %do;
        %put ERROR: result_export: failed while writing OC.;
        %return;
    %end;


    %put NOTE: result_export completed successfully.;
    %put NOTE: result_export workbook=&_outfile.;

%mend result_export;