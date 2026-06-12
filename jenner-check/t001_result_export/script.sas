/* Exercises the %result_export utility from sas/blrm_export.sas, which writes
   the BLRM simulation outputs to an Excel workbook. Reproduced as the author
   wrote it; the only change is the path separator in the OUTFILE= (Windows
   "\" -> portable "/") so the workbook lands in the working directory. The
   six PROC EXPORT steps, DBMS=XLSX, the SHEET= names and the dataset names
   are unchanged. */

%macro result_export(path=,file_name=,info=);
proc export data=&info.
	outfile="&path./&file_name..xlsx"
	dbms=xlsx
	replace;
	sheet= "INFO";
run;



proc export data=dlt_all
	outfile="&path./&file_name..xlsx"
	dbms=xlsx
	replace;
	sheet= "DLT";
run;
proc export data=EWOC_all
	outfile="&path./&file_name..xlsx"
	dbms=xlsx
	replace;
	sheet= "EWOC";
run;

proc export data=SIM_all
	outfile="&path./&file_name..xlsx"
	dbms=xlsx
	replace;
	sheet= "SIM";
run;
proc export data=DOSE_SELECTION
	outfile="&path./&file_name..xlsx"
	dbms=xlsx
	replace;
	sheet= "DOSE_SELECTION";
run;
proc export data=SIM_RESULTS
	outfile="&path./&file_name..xlsx"
	dbms=xlsx
	replace;
	sheet= "SIM_RESULTS";
run;
%mend;

%result_export(path=., file_name=blrm_results, info=info);
