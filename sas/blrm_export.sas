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
%mend;

%macro result_conti_export(path=,file_name=,info=);
data &info.;
	set &info.;
	ewoc=&ewoc.;
run;
proc export data=&info.
	outfile="&path.\&file_name..xlsx"
	dbms=xlsx
	replace;
	sheet= "INFO";
run;
	


proc export data=dlt_record_all
	outfile="&path.\&file_name..xlsx"
	dbms=xlsx
	replace;
	sheet= "DLT";
run;
proc export data=EWOC_record_all
	outfile="&path.\&file_name..xlsx"
	dbms=xlsx
	replace;
	sheet= "EWOC";
run;

proc export data=SIM_record_all
	outfile="&path.\&file_name..xlsx"
	dbms=xlsx
	replace;
	sheet= "SIM";
run;

%mend;