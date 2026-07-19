/*======================================================================*/
/* result_export                                                        */
/*                                                                      */
/* path      : output directory                                         */
/* file_name : output workbook name without the .xlsx extension         */
/* info      : analysis-level information data set                      */
/*                                                                      */
/* The macro also exports the standard data sets produced by the        */
/* current discrete-dose BLRM analysis and simulation workflow.         */
/*======================================================================*/
%macro result_export(path=,file_name=,info=);
proc export data=&info.
	outfile="&path.\&file_name..xlsx"
	dbms=xlsx
	replace;
	sheet= "INFO";
run;

proc export data=dlt_all
	outfile="&path.\&file_name..xlsx"
	dbms=xlsx
	replace;
	sheet= "DLT";
run;

proc export data=EWOC_all
	outfile="&path.\&file_name..xlsx"
	dbms=xlsx
	replace;
	sheet= "EWOC";
run;

proc export data=SIM_all
	outfile="&path.\&file_name..xlsx"
	dbms=xlsx
	replace;
	sheet= "SIM";
run;

proc export data=DOSE_SELECTION
	outfile="&path.\&file_name..xlsx"
	dbms=xlsx
	replace;
	sheet= "DOSE_SELECTION";
run;

proc export data=SIM_RESULTS
	outfile="&path.\&file_name..xlsx"
	dbms=xlsx
	replace;
	sheet= "SIM_RESULTS";
run;
%mend result_export;
