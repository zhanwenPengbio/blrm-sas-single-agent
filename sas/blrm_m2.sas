/*blrm_main: the main part*/

/*datain: the input dataset.*/
/*dataout: the name of output dataset.*/
/*ref_dose: reference dose, for scaling the dose range.*/
/*mu1: the prior mean of log(alpha)*/
/*mu2: the prior mean of log(beta)*/
/*v1: the prior variance of log(alpha)*/
/*v2: the prior variance of log(beta)*/
/*rho: the prior covariance of log(alpha) and log(beta)*/
/*simming: simulation indicator, if in simulation, postSumint will not be output*/
%macro blrm_main(datain=,dataout=, ref_dose=, mu1=, mu2=, v1=, v2=, rho= ,
alpha=0.05, nbi=2000, nmc = 20000, simming=0,seed = 9527);
/*ods graphics off;*/
ods exclude all;
%if &simming = 0 %then %do;
ods exclude close;
ods output PostSumInt=PostSumInt;
%end;
proc mcmc data=&datain. stats(alpha=&alpha.) seed=&seed. nbi=&nbi. nmc=&nmc.
	outpost=&dataout. ;

	array log_par[2] log_alpha log_beta;
	array mu_prior[2];
	array Sigma[2,2];

	begincnst;
		call zeromatrix(mu_prior);
		call identity(Sigma);
		mu_prior[1] = &mu1.;
		mu_prior[2] = &mu2.;
		Sigma[1,1] = &v1.;
		Sigma[2,2] = &v2.;
		Sigma[1,2] = &rho.;
		Sigma[2,1] = &rho.;
	endcnst;
	
	prior log_par ~ mvn(mu_prior,Sigma);
	
	parms log_par;
	
	delta = log_par[1] + exp(log_par[2])*log(dose/&ref_dose.);

	pi_d = logistic(delta);

	model dltn ~ binomial(n,pi_d);

/*	preddist outpred=pout  nsim=&nsim.;*/
run;
ods graphics off;
%mend blrm_main;

%macro loop(values);    
                                                                                                                
     /* Count the number of values in the string */                                                                                                                                   
     %let count=%sysfunc(countw(&values)); 

     /* Loop through the total number of values */                                                                                         
     %do i = 1 %to &count;                                                                                                              
      %let value=%qscan(&values,&i,%str(,));                                                                                            
      %put &value;                                                                                                                      
     %end;                                                                                                                              
                                                                                                                                        
%mend;  



/*blrm_stat: a macro that summarize posterior dlt from blrm_main.*/
/*dose_list: the dose levels that is wanted to summarize the posterior dlt, not necessary the same dose levels in blrm_main.*/
/*(example of dose_list: %str(2,4,5,12))*/
/*mcmcout: the output dataset from*/
/*dataout: the summary of blrm output*/
/*ref_dose: reference dose, for scaling the dose range.*/
/*stat: the statistics want to summarize*/
/*debug: if set to 0, temporary dataset will be deleted.*/

%macro blrm_stat(dose_list=%str(),ref_dose=,mcmcout=,stat=%str(mean, std, q1, median, q3),dataout=,debug=0);

%let count=%sysfunc(countw(&dose_list.)); 
data temp1;
	set &mcmcout.;
	%do i = 1 %to &count.;
		%let vt=%qscan(&dose_list.,&i,%str(,));
		%LET vx = %sysfunc(cats(dose,_,&vt.));
		%put &=vt.;
		%put &=vx.;
		&vx.=logistic(log_alpha + exp(log_beta)*log(&vt./&ref_dose.));
	%end;
run;
%let count2=%sysfunc(countw(&stat.)); 
data stat_temp;
	%do j=1 %to &count2.;
		%let vt=%qscan(&stat.,&j,%str(,));
	g1 = cats("&vt.",'=');output;
	%end;
run;
proc sql noprint;
	select g1 into: stats separated by ' '
	from stat_temp
;quit;
proc means data= temp1 noprint;
	var dose:;
	output out=temp2(drop=_type_ _freq_) &stats./ autoname; 

run;

proc transpose data=temp2 out=temp3;
run;

data temp4;
	set temp3;
	varname=cats(scan(_name_,1,'_'),':',scan(_name_,2,'_'));
	stat=scan(_name_,3,'_');
	drop _name_;
run;

proc sort data=temp4;
	by varname;
run;

proc transpose data=temp4 out=temp5(drop=_name_ );
	by varname;
	id stat;
	var col1;
run;

data &dataout.;
	retain dose ;
	set temp5;
	dose=input(compress(varname,':','a'),best.);
	drop varname;
	proc sort;
		by dose;
run;

%if &debug.=0 %then %do;
proc delete data=temp1-temp5 stat_temp;
%end;

%mend blrm_stat;

/*blrm_ewoc: a macro that summarize under-dose, target dose, over dose.*/
/*dose_list: the dose levels that is wanted to summarize the posterior dlt, not necessary the same dose levels in blrm_main.*/
/*(example of dose_list: %str(2,4,5,12))*/
/*mcmcout: the output dataset from*/
/*dataout: the summary of blrm output*/
/*ud: under dosing boundary*/
/*od: over dosing boundary*/
/*ewoc: the escalation with overdose control RULE, control the posterior dlt not exceed some specific value*/
/*debug: if set to 0, temporary dataset will be deleted.*/

%macro blrm_ewoc(dose_list=,ref_dose=,mcmcout=out_mcmc,dataout=d1,ud=0.16,od=0.33,ewoc=0.25,debug=0);
%let count=%sysfunc(countw(&dose_list.)); 
data etemp1;
	set &mcmcout.;
	%do i = 1 %to &count.;
		%let vt=%qscan(&dose_list.,&i,%str(,));
		%LET vx = %sysfunc(cats(dose,_,&vt.));
		&vx.=logistic(log_alpha + exp(log_beta)*log(&vt./&ref_dose.));;
	%end;
run;

data etemp2;
	set etemp1;
	%do i = 1 %to &count.;
		%let vt=%qscan(&dose_list.,&i,%str(,));
		%LET vx = %sysfunc(cats(dose,_,&vt.));
/*		%put &=vx.;*/
		if &vx.<= &ud. then &vx._ud = 1;
		else if &vx.> &ud. and &vx.<= &od. then &vx._tt = 1;
		else if &vx.> &od. then &vx._od = 1;
	%end;
run;
proc contents data=etemp2 out=var_list(keep=name) noprint ;run;
proc sql noprint ;
	select name into :udvars separated by ' '
from var_list where upcase(name) like '%^_UD' escape '^';
	select name into :ttvars separated by ' '
from var_list where upcase(name) like '%^_TT' escape '^';
	select name into :odvars separated by ' '
from var_list where upcase(name) like '%^_OD' escape '^';
;
quit;
proc delete data=var_list;
data etemp3;
	set etemp2;
	array xx(*) &udvars. &ttvars. &odvars.;
	do i = 1 to dim(xx);
		if missing(xx(i)) then xx(i) = 0;
	end;
	drop i;
run;

/*%put &=udvars.;*/
/*%put &=ttvars.;*/
/*%put &=odvars.;*/
Proc means data=etemp3  maxdec=0 noprint;

   var &udvars. &ttvars. &odvars.;
output out=etemp4(drop=_:) mean=/autoname;
run;

proc transpose data=etemp4 out=etemp5;
run;

data etemp6;
	set etemp5;
	varname=cats(scan(_name_,1,'_'),':',scan(_name_,2,'_'));
	/*stat='Prob';*/
	stat=catx('_',upcase(scan(_name_,3,'_')),scan(_name_,4,'_'));
	drop _name_;
run;

proc sort data=etemp6;
	by varname ;
run;

proc transpose data=etemp6 out=etemp7(drop=_name_ );
	by varname;
	id stat;
	var COL1;
run;

data &dataout.;
	retain dose;
	set etemp7;
	dose=input(compress(varname,':','a'),best.);
	if od_mean > &ewoc. then EWOC = 'False';
	else EWOC = 'True';
	%if &debug.=1 %then %do;
		check1 = sum(ud_mean,tt_mean,od_mean);
	%end;
	drop varname;
	proc sort;
		by dose;
/*	rename varname = dose_info;*/
run;




%if &debug.=0 %then %do;
	proc delete data=etemp1-etemp7;
%end;

%mend blrm_ewoc;

/*blrm_sim_one: a macro that perform 1 simulation of a blrm trial.*/
/*use macros: blrm_main, blrm_ewoc*/
/*stop rule:*/
/*1. All sample size are used.*/
/*2. All dose levels not satisified EWOC.*/
/*other params:*/
/*true_tox: dataset that contains the true toxicity of each dose level.*/
/*start_dose: the starting dose of simulation, may not necessary be the lowest dose level*/
/*max_sample: the maximum sample size*/
/*cohortn: the cohort size of each enrollment, if max_sample mod(cohortn) ^= 0, then the remainder will be the size of the last cohort.*/
/*skip: Default is 0, which means can only select the next dose level when escalating. If set to 1, skipping dose is allowed. ;*/
%macro blrm_sim_one(dose_list=%str(),ref_dose=,mu1=, mu2=, v1=, v2=, rho= ,
        alpha=0.05, nbi=1000, nmc = 10000, 
        ud=0.16,od=0.33,ewoc=0.25,
        true_tox=,start_dose=,max_sample=,cohortn=,skip=,seed=,sim_batch=,debug=0);

    /* PREFIXED LOCAL VARIABLES */
    %local bs1_ebn bs1_iter_count bs1_max_iter 
           bs1_cur_dose bs1_cur_ord bs1_cur_tox bs1_ewoc_left;
    
    %let bs1_ebn = 1;
    %let bs1_sample_expand = 0;
    %let bs1_max_iter = %eval(&max_sample./&cohortn. + 20);
    %let bs1_iter_count = 1;

    proc sort data=&true_tox.;
        by dose;
    run;

    data dose_sort;
        set &true_tox.;
        by dose;
        ord = _n_;
        keep dose ord;
    run;

    data dlt_record;
        set _null_;
    run;
    
    data ewoc_record;
        set _null_;
    run;
    
    %let bs1_cur_dose = &start_dose.;
    proc sql noprint;
        select ord into: bs1_cur_ord trimmed from dose_sort
        where dose = &bs1_cur_dose.;
    quit;
    
    %do %until (&bs1_sample_expand>=&max_sample. or &bs1_iter_count. > &bs1_max_iter.);
        /* Cohort size adjustment */
        %if %eval(&max_sample.-&bs1_sample_expand.)<=&cohortn. and %eval(&max_sample.-&bs1_sample_expand.)>0 %then %do;
            %let cohortn = %eval(&max_sample.-&bs1_sample_expand.);
        %end;
        
        /* Get current toxicity - with null check */
        %if %sysevalf(%superq(bs1_cur_dose)=,boolean) %then %goto bs1_exit;
        proc sql noprint;
            select true_tox into: bs1_cur_tox trimmed
            from &true_tox. where dose = &bs1_cur_dose.;
        quit;
        
        /* Generate cohort data */
        data tox_gen;
            dose = &bs1_cur_dose.;
            n = &cohortn.;
            call streaminit(%eval(&seed. + &sim_batch.*10000 + &bs1_ebn.*100)); 
            dltn = rand("BINOMIAL", &bs1_cur_tox., n);
            true_tox = &bs1_cur_tox.;
            enroll_batch = &bs1_ebn.;
        run;
        
        /* Update records */
        data dlt_record;
            set dlt_record tox_gen;
        run;
        
        /* Process BLRM */
        proc sql;
            create table tox_sum as
            select dose, sum(n) as n, sum(dltn) as dltn
            from dlt_record group by dose;
        quit;
        
        %blrm_main(datain=tox_sum, dataout=out_mcmc, ref_dose=&ref_dose., 
            mu1=&mu1., mu2=&mu2., v1=&v1., v2=&v2., rho=&rho.,
            alpha=&alpha., nbi=&nbi., nmc=&nmc., 
            seed=%eval(&seed.+&bs1_ebn.), simming=1);
        
        %blrm_ewoc(dose_list=&dose_list., ref_dose=&ref_dose., mcmcout=out_mcmc,
            dataout=next_dose, ud=&ud., od=&od., ewoc=&ewoc., debug=&debug.);
        
        /* Update EWOC records */
        data ewoc_temp;
            set next_dose;
            enroll_batch = &bs1_ebn.;
        run;
        
        data ewoc_record;
            set ewoc_record ewoc_temp;  
        run;

        /* Dose selection logic */
        data tt1;
            merge dose_sort next_dose;
            by dose;
            if EWOC = 'True';
            %if &skip=0 %then %do;
                if ord <= &bs1_cur_ord.+1;
            %end;
            %else %if &skip=1 %then %do;
                if ord <= &bs1_cur_ord.+2;
            %end;
            proc sort;
                by descending TT_Mean;
        run;
        
        /* EWOC validation */
        data _null_;
            if 0 then set tt1 nobs=nobs;
            call symputx('bs1_ewoc_left', nobs);
            stop;
        run;
        
		data tt2;
			set _null_;run;


        %if &bs1_ewoc_left = 0 %then %do;
			%let bs1_cur_dose = .;
			%let bs1_sample_expand = %eval(&bs1_sample_expand. + &cohortn.);
			%goto bs1_exit;
		%end;
        
        /* Select next dose */
        data tt2;
            set tt1;
            by descending TT_Mean;
            if _n_ = 1;
            if not missing(dose) then do;
                call symputx('bs1_cur_dose', dose);
                call symputx('bs1_cur_ord', ord);
            end;
        run;
        
        /* Update counters */
        %let bs1_sample_expand = %eval(&bs1_sample_expand. + &cohortn.);
        %let bs1_ebn = %eval(&bs1_ebn.+1);
        %let bs1_iter_count = %eval(&bs1_iter_count.+1);
    %end;

%bs1_exit:
    /* Final output processing */
    data dlt_record;
        set dlt_record;
        sim_batch = &sim_batch.;
    run;
    
    data ewoc_record;
        set ewoc_record;
        sim_batch = &sim_batch.;
    run;
    
    data sim_record;
        dose_select = &bs1_cur_dose.;
        sample_used = &bs1_sample_expand.;
        sim_batch = &sim_batch.;
    run;
    
    %if &debug. = 0 %then %do;
        proc delete data=DOSE_SORT TOX_GEN TOX_SUM OUT_MCMC NEXT_DOSE EWOC_TEMP tt1 tt2;
    %end;
%mend blrm_sim_one;

%macro blrm_sim_n(dose_list=%str(),ref_dose=,true_tox=,start_dose=,max_sample=,cohortn=,skip=,
        mu1=, mu2=, v1=, v2=, rho= ,
        alpha=0.05, nbi=1000, nmc = 10000, 
        ud=0.16,od=0.33,ewoc=0.25,
        seed=,sim_time=4,debug=0, time_out=time_summary,time_sum=1);
        
    %local bs_n_i start_time end_time total_time;
    data dlt_all;
        set _null_;
    run;
    data ewoc_all;
        set _null_;
    run;
    data sim_all;
        set _null_;
    run;
    data &time_out;
        set _null_;
    run;
    
    /* Initialize timing dataset */
    data time_record;
        sim_batch = .;
        start_dt = .;
        end_dt = .;
        elapsed_sec = .;
        delete;
    run;
    
    %let total_start = %sysfunc(datetime());
    
    %do bs_n_i=1 %to &sim_time.;
        %let sim_start = %sysfunc(datetime());
        
        %blrm_sim_one(dose_list=&dose_list.,
                ref_dose=&ref_dose.,mu1=&mu1., mu2=&mu2., v1=&v1., v2=&v2., rho=&rho. ,
                alpha=&alpha., nbi=&nbi., nmc = &nmc.,
                ud=&ud.,od=&od.,ewoc=&ewoc.,
                true_tox=&true_tox.,
                start_dose=&start_dose.,
                max_sample=&max_sample,cohortn=&cohortn.,skip=&skip.,
                seed=%eval(&seed.+1024+&bs_n_i.*9),
                sim_batch=&bs_n_i.,
                debug=&debug.);
        	
        %let sim_end = %sysfunc(datetime());
        
        /* Record individual simulation time */
        data temp_time;
            sim_batch = &bs_n_i;
            start_dt = &sim_start;
            end_dt = &sim_end;
            elapsed_sec = &sim_end - &sim_start;
/*            output;*/
        run;
        
        /* Append to time record */
        proc append base=time_record data=temp_time;
        run;
        
        data dlt_all;
            set dlt_all dlt_record;
        run;
        data ewoc_all;
            set ewoc_all ewoc_record;
        run;
        data sim_all;
            set sim_all sim_record;
        run;
		%if &debug=0 %then %do;
			proc delete data=dlt_record ewoc_record sim_record;
			%end;
    %end;
    
    %let total_end = %sysfunc(datetime());
    
    /* Calculate total time */
    data &time_out;
        total_start = &total_start;
        total_end = &total_end;
        total_elapsed_sec = total_end - total_start;
        avg_sec_per_sim = total_elapsed_sec / &sim_time.;
        format total_start total_end datetime19.;
    run;
    
    /* Append individual times to summary */
    data time_record;
        set time_record;
        format start_dt end_dt datetime19.;
    run;
    
    /* Combine with total time */
    data &time_out;
        set &time_out time_record;
    run;
    
    /* Print timing summary */
	%if &time_sum=1 %then %do;
	ods exclude close;
    title "Simulation Timing Summary";
/*    proc print data=&time_out;*/
/*        var sim_batch start_dt end_dt elapsed_sec;*/
/*        where not missing(sim_batch);*/
/*    run;*/
    proc print data=&time_out;
        where missing(sim_batch);
        var total_start total_end total_elapsed_sec avg_sec_per_sim;
    run;
    title;
	%end;
%mend blrm_sim_n;

%macro blrm_sim_sum(
    sim_data=,              /* Output from blrm_sim_n (sim_all) */
    dlt_data=,              /* Patient-level data from blrm_sim_n (dlt_all) */
    true_tox=,              /* Dataset with true toxicity probabilities */
    target_range_low=0.16,  /* Lower bound of target toxicity range */
    target_range_high=0.33, /* Upper bound of target toxicity range (also overdose threshold) */
    acceptable_low=0.05,    /* Lower bound for acceptable dose */
    out_summary=  ,          /* Output dataset name */
	debug=0
);

/* Use target_range_high as overdose threshold */
%let over_tox_threshold = &target_range_high;

/* Step 1: Identify doses in different categories */
proc sql;
    create table dose_categories as
    select 
        dose,
        true_tox,
        /* Target dose: within target range */
        case when true_tox between &target_range_low and &target_range_high 
             then 1 else 0 end as is_target,
        /* Acceptable dose: above min acceptable but not over toxic */
        case when true_tox between &acceptable_low and &target_range_high 
             then 1 else 0 end as is_acceptable,
        /* Over toxic dose: above target range high */
        case when true_tox > &target_range_high 
             then 1 else 0 end as is_over_toxic
    from &true_tox;
quit;

/* Step 2: Calculate dose selection proportions */
proc sql;
    create table dose_selection as
    select 
        a.sim_batch,
        a.dose_select,
        max(b.is_target) as is_target,
        max(b.is_acceptable) as is_acceptable,
        max(b.is_over_toxic) as is_over_toxic
    from &sim_data a
    left join dose_categories b on a.dose_select = b.dose
    group by a.sim_batch, a.dose_select;
quit;

proc means data=dose_selection noprint nway;
    var is_target is_acceptable is_over_toxic;
    output out=selection_summary mean= / autoname;
run;

/* Step 3: Calculate patients exposed to over-toxic doses */
proc sql;
    create table over_toxic_exposure as
    select 
        d.sim_batch,
        sum(case when t.true_tox > &target_range_high then d.n else 0 end) as n_over_toxic
    from &dlt_data d
    left join &true_tox t on d.dose = t.dose
    group by d.sim_batch;
quit;

proc means data=over_toxic_exposure noprint;
    var n_over_toxic;
    output out=exposure_summary mean=avg_n_over_toxic;
run;

/* Step 4: Combine all summary metrics */
data &out_summary;
    length Metric $40 Value 8;
    
    /* Selection proportions */
    set selection_summary(keep=is_target_mean 
                             is_acceptable_mean 
                             is_over_toxic_mean);
    
    Metric = "Prop_Correct_Dose"; 
    Value = is_target_mean; 
    output;
    
    Metric = "Prop_Acceptable_Dose"; 
    Value = is_acceptable_mean; 
    output;
    
    Metric = "Prop_Over_Toxic_Dose"; 
    Value = is_over_toxic_mean; 
    output;
    
    /* Exposure metric */
    set exposure_summary(keep=avg_n_over_toxic);
    Metric = "Avg_Patients_Over_Toxic";
    Value = avg_n_over_toxic;
    output;
    
    keep Metric Value;
    
    /* Add parameters to output */
    call symput('accept_low', &acceptable_low);
    call symput('target_low', &target_range_low);
    call symput('target_high', &target_range_high);
run;

/* Step 5: Add parameter information to dataset */
data &out_summary.;
    set &out_summary;
	format _all_;
	length Parameters $170.;
    Parameters = cats("Acceptable: >=", &acceptable_low., " | Target: ", &target_range_low., 
                      "-", &target_range_high., " | Over: >", &target_range_high.);
/*    label */
/*        Parameters = "Dose Category Definitions"*/
/*        Metric = "Performance Metric"*/
/*        Value = "Result Value";*/
run;

/* Step 6: Clean up */
 %if &debug=0 %then %do;
proc datasets library=work nolist;
    delete dose_categories dose_selection selection_summary
           over_toxic_exposure exposure_summary ;
quit;
%end;
%mend blrm_sim_sum;

%macro blrm_sim_conti_one(
    /* Design & BLRM */
    dose_list=%str(5,10,20,40,80),
    ref_dose=,
    mu1=, mu2=, v1=, v2=, rho=,
    alpha=0.05, nbi=1000, nmc=10000,
    ud=0.16, od=0.33, ewoc=0.25,

    /* Trial logistics */
    start_dose=,
    max_sample=,
    cohortn=3,
    skip=0,
    seed=12345,
    sim_batch=1,

    /* Truth via FCMP */
    tox_fn=,
    scenario=Medium,

    /* Uniform disturbance on ACTUAL dose (per patient) */
    disturb_actual_lo=-1, disturb_actual_hi=1,

    /* Out-of-range p handling */
    abort_on_bad_p=0,

    /* Misc */
    debug=0
  );
  %if %sysevalf(&disturb_actual_hi < &disturb_actual_lo) %then %do;
    %put ERROR: disturb_actual_hi < disturb_actual_lo.;
    %return;
  %end;
    /* PREFIXED LOCAL VARIABLES */
    %local bs1_ebn bs1_iter_count bs1_max_iter 
           bs1_cur_dose bs1_cur_ord bs1_cur_tox bs1_ewoc_left;
    
    %let bs1_ebn = 1;
    %let bs1_sample_expand = 0;
    %let bs1_max_iter = %eval(&max_sample./&cohortn. + 20);
    %let bs1_iter_count = 1;
	%let doses_sp = %sysfunc(compbl(%sysfunc(tranwrd(&dose_list., %str(,), %str( )))));
	%let n_dose   = %sysfunc(countw(&doses_sp));
	%let step_limit = %sysfunc(ifc(&skip,2,1));

/*    proc sort data=&true_tox.;*/
/*        by dose;*/
/*    run;*/

    data dose_sort;
	    length scenario $40;
	    scenario="&scenario";
	    %do i=1 %to &n_dose;
	      ord=&i;
	      dose = input(scan("&doses_sp", &i, ' '), best.);
	      output;
	    %end;
    run;

	
    data dlt_record; length scenario $40;  stop; run;
    data ewoc_record;  length scenario $40; stop; run;

/*	since we use distrubance and continuous toxicity curve, p>1 situation could happen*/
    data bad_p_events; length scenario $40 msg $200; stop; run;

    %let bs1_cur_dose = &start_dose.;
    proc sql noprint;
        select ord into: bs1_cur_ord trimmed from dose_sort
        where dose = &bs1_cur_dose.;
    quit;
    
    %do %until (&bs1_sample_expand>=&max_sample. or &bs1_iter_count. > &bs1_max_iter.);
        /* Cohort size adjustment */
        %if %eval(&max_sample.-&bs1_sample_expand.)<=&cohortn. and %eval(&max_sample.-&bs1_sample_expand.)>0 %then %do;
            %let cohortn = %eval(&max_sample.-&bs1_sample_expand.);
        %end;
        
        
        /* Generate cohort data */
/*		Even in same dose level, different enrolled subject may have different actual dose*/
/*		leading to different actual toxicity in the curve.*/
		data tox_gen;
			scenario="&scenario";
			plan_dose = &bs1_cur_dose.;
			call streaminit(%eval(&seed + 200000*&sim_batch. + 10000*&bs1_ebn.));
			do i=1 to &cohortn.;
				pid=i;
				jitter = &disturb_actual_lo + (&disturb_actual_hi - &disturb_actual_lo)*rand('uniform');
				act_dose = plan_dose + jitter;
				%if %sysevalf(%superq(tox_fn)=,boolean) %then %do;
          			p_raw = p_true(act_dose, "&scenario");
					plan_tox = p_true(plan_dose, "&scenario");
        		%end;
        		%else %do;
          			p_raw = &tox_fn.(act_dose);
					plan_tox = &tox_fn.(plan_dose);
       			%end;

			if missing(p_raw) or p_raw<0 or p_raw>1 then do;
	          bad_count + 1;
	          msg = catx(' ',
	                'Toxicity function returned invalid p:',
	                put(p_raw, best16.),
	                'at actual=', put(act_dose, best16.),
	                'planned=', put(plan_dose, best16.),
	                '(sim_batch=&sim_batch, batch=&bs1_ebn., pid=', put(pid, 8.), ')');
	          putlog 'ERROR: ' msg;
	          if &abort_on_bad_p then call symputx('_abort_flag', 1);
	        end;
			act_tox = p_raw; if act_tox<0 then act_tox=0; else if act_tox>1 then act_tox=1;
			y = rand('bernoulli', act_tox);
			enroll_batch = &bs1_ebn.;
/*			delete i;*/
			output;
			end;
	    run;

        
        /* Update records */
        data dlt_record;
            set dlt_record tox_gen;
        run;
        
        /* Process BLRM */
/*		in PROC MCMC, the input must be all the DLT record, NOT each enrollment*/
        proc sql;
            create table tox_sum as
            select plan_dose as dose, count(*) as n, sum(y) as dltn
            from dlt_record group by plan_dose;
        quit;
        
        %blrm_main(datain=tox_sum, dataout=out_mcmc, ref_dose=&ref_dose., 
            mu1=&mu1., mu2=&mu2., v1=&v1., v2=&v2., rho=&rho.,
            alpha=&alpha., nbi=&nbi., nmc=&nmc., 
            seed=%eval(&seed.+&bs1_ebn.), simming=1);
        
        %blrm_ewoc(dose_list=&dose_list., ref_dose=&ref_dose., mcmcout=out_mcmc,
            dataout=next_dose, ud=&ud., od=&od., ewoc=&ewoc., debug=&debug.);
        
        /* Update EWOC records */
        data ewoc_temp;
            set next_dose;
            enroll_batch = &bs1_ebn.;
        run;
        
        data ewoc_record;
            set ewoc_record ewoc_temp;  
        run;

        /* Dose selection logic */
        data tt1;
            merge dose_sort next_dose;
            by dose;
            if EWOC = 'True';
            %if &skip=0 %then %do;
                if ord <= &bs1_cur_ord.+1;
            %end;
            %else %if &skip=1 %then %do;
                if ord <= &bs1_cur_ord.+2;
            %end;
            proc sort;
                by descending TT_Mean;
        run;
        
        /* EWOC validation */
        data _null_;
            if 0 then set tt1 nobs=nobs;
            call symputx('bs1_ewoc_left', nobs);
            stop;
        run;
        
		data tt2;
			set _null_;run;


        %if &bs1_ewoc_left = 0 %then %do;
			%let bs1_cur_dose = .;
			%let bs1_sample_expand = %eval(&bs1_sample_expand. + &cohortn.);
			%goto bs1_exit;
		%end;
        
        /* Select next dose */
        data tt2;
            set tt1;
            by descending TT_Mean;
            if _n_ = 1;
            if not missing(dose) then do;
                call symputx('bs1_cur_dose', dose);
                call symputx('bs1_cur_ord', ord);
            end;
        run;
        
        /* Update counters */
        %let bs1_sample_expand = %eval(&bs1_sample_expand. + &cohortn.);
        %let bs1_ebn = %eval(&bs1_ebn.+1);
        %let bs1_iter_count = %eval(&bs1_iter_count.+1);
    %end;

%bs1_exit:
    /* Final output processing */
    data dlt_record;
        set dlt_record;
        sim_batch = &sim_batch.;
    run;
    
    data ewoc_record;
        set ewoc_record;
        sim_batch = &sim_batch.;
    run;
    
    data sim_record;
        dose_select = &bs1_cur_dose.;
        sample_used = &bs1_sample_expand.;
        sim_batch = &sim_batch.;
    run;
    
    %if &debug. = 0 %then %do;
        proc delete data=DOSE_SORT TOX_GEN TOX_SUM OUT_MCMC NEXT_DOSE EWOC_TEMP tt1 tt2;
    %end;

%mend blrm_sim_conti_one;

/*==============================================================
/*==============================================================
  Batch runner with timing summary for continuous-dose BLRM
  Depends on: %blrm_sim_conti_one (your inner macro)
  Produces:
    dlt_record_all, ewoc_record_all, sim_record_all, bad_p_events_all
    sim_summary, sim_selection_dist
    <time_out> (default: time_summary) with total + per-sim timings
  ==============================================================*/
/*==============================================================
  BLRM continuous-dose simulation (per-sim) + batch runner w/ timing
  - %blrm_sim_conti_one : one simulation, per-patient jitter (actual dose)
  - %blrm_sim_conti_n   : runs N sims, aggregates outputs, writes timing
  NOTE: Requires your %blrm_main and %blrm_ewoc to be defined.
  ==============================================================*/

/*-------------------------------*
 | One simulation (continuous)   |
 *-------------------------------*/
%macro blrm_sim_conti_one(
    /* Design & BLRM */
    dose_list=%str(5,10,20,40,80),
    ref_dose=,
    mu1=, mu2=, v1=, v2=, rho=,
    alpha=0.05, nbi=1000, nmc=10000,
    ud=0.16, od=0.33, ewoc=0.25,

    /* Trial logistics */
    start_dose=,
    max_sample=,
    cohortn=3,
    skip=0,
    seed=12345,
    sim_batch=1,

    /* Truth via FCMP (or p_true if tox_fn blank) */
    tox_fn=,
    scenario=Medium,

    /* Uniform disturbance on ACTUAL dose (per patient) */
    disturb_actual_lo=-1, disturb_actual_hi=1,

    /* Out-of-range p handling */
    abort_on_bad_p=0,

    /* Misc */
    debug=0
  );

  /* input checks */
  %if %sysevalf(&disturb_actual_hi < &disturb_actual_lo) %then %do;
    %put ERROR: disturb_actual_hi < disturb_actual_lo.;
    %return;
  %end;

  /* Local state */
  %local bs1_ebn bs1_iter_count bs1_max_iter bs1_sample_expand
         bs1_cur_dose bs1_cur_ord bs1_cur_tox bs1_ewoc_left
         doses_sp n_dose step_limit _abort_flag;

  %let bs1_ebn = 1;
  %let bs1_sample_expand = 0;
  /* allow headroom beyond max_sample */
  %let bs1_max_iter = %sysfunc(ceil(%sysevalf(&max_sample./&cohortn.))) + 20;
  %let bs1_iter_count = 1;
  %let doses_sp = %sysfunc(compbl(%sysfunc(tranwrd(&dose_list., %str(,), %str( )))));
  %let n_dose   = %sysfunc(countw(&doses_sp));
  %let step_limit = %sysfunc(ifc(&skip,2,1));
  %let _abort_flag=0;

  /* Design dose order */
  data dose_sort;
    length scenario $40;
    scenario="&scenario";
    %do i=1 %to &n_dose;
      ord=&i;
      dose = input(scan("&doses_sp", &i, ' '), best.);
      output;
    %end;
  run;

  /* Output shells */
  data dlt_record;    length scenario $40; stop; run;
  data ewoc_record;   length scenario $40; stop; run;
  data bad_p_events;  length scenario $40 msg $200; stop; run;

  /* Starting dose -> ord (guard exact match) */
  %let bs1_cur_dose = &start_dose.;
  proc sql noprint;
    select ord into :bs1_cur_ord trimmed
    from dose_sort
    where dose = &bs1_cur_dose.;
  quit;
  %if %sysevalf(%superq(bs1_cur_ord)=,boolean) %then %do;
    %put ERROR: start_dose=&start_dose not found in dose_list=&dose_list.;
    %return;
  %end;

  /*==================== enrollment loop ====================*/
  %do %until (&bs1_sample_expand>=&max_sample. or &bs1_iter_count. > &bs1_max_iter.);

    /* last cohort size adjustment */
    %if %eval(&max_sample.-&bs1_sample_expand.)<=&cohortn. and %eval(&max_sample.-&bs1_sample_expand.)>0 %then %do;
      %let cohortn = %eval(&max_sample.-&bs1_sample_expand.);
    %end;

    /* Generate cohort (patient-level jitter on actual dose) */
    data tox_gen;
      length scenario $40 msg $200;
      scenario="&scenario";
      plan_dose = &bs1_cur_dose.;
      call streaminit(%eval(&seed + 200000*&sim_batch. + 10000*&bs1_ebn.));
      do i=1 to &cohortn.;
        pid=i;
        jitter   = &disturb_actual_lo + (&disturb_actual_hi - &disturb_actual_lo)*rand('uniform');
        act_dose = plan_dose + jitter;

        /* If your tox_fn/p_true expects in-range values only, optionally clip:
           act_dose = max(min(act_dose, max(&doses_sp)), min(&doses_sp)); */

        %if %sysevalf(%superq(tox_fn)=,boolean) %then %do;
          p_raw   = p_true(act_dose, "&scenario");
          plan_tox= p_true(plan_dose, "&scenario");
        %end;
        %else %do;
          p_raw   = &tox_fn.(act_dose);
          plan_tox= &tox_fn.(plan_dose);
        %end;

        if missing(p_raw) or p_raw<0 or p_raw>1 then do;
          bad_count + 1;
          msg = catx(' ',
                 'Invalid p:', put(p_raw, best16.),
                 'at actual=', put(act_dose, best16.),
                 'planned=', put(plan_dose, best16.),
                 '(sim_batch=&sim_batch, batch=&bs1_ebn., pid=', put(pid, 8.), ')');
          putlog 'ERROR: ' msg;
          scenario="&scenario";
          output; /* keep this row in tox_gen even if invalid; will be clipped below */
          if &abort_on_bad_p then call symputx('_abort_flag', 1);
        end;
        else output;

      end;
      drop i msg;
    run;

    %if %sysevalf(&_abort_flag) %then %do;
      %put ERROR: Aborting this simulation due to invalid toxicity p (abort_on_bad_p=1).;
      %goto bs1_exit;
    %end;

    /* clamp probabilities to [0,1] and realize DLT */
    data tox_gen;
      set tox_gen;
      act_tox = max(0, min(1, p_raw));
      y = rand('bernoulli', act_tox);
      enroll_batch = &bs1_ebn.;
    run;

    /* Append patient records */
    data dlt_record; set dlt_record tox_gen; run;

    /* Summarize for BLRM (by plan_dose) */
    proc sql;
      create table tox_sum as
      select plan_dose as dose, count(*) as n, sum(y) as dltn
      from dlt_record
      group by plan_dose;
    quit;

    /* Run BLRM and EWOC */
    %blrm_main(datain=tox_sum, dataout=out_mcmc, ref_dose=&ref_dose., 
      mu1=&mu1., mu2=&mu2., v1=&v1., v2=&v2., rho=&rho.,
      alpha=&alpha., nbi=&nbi., nmc=&nmc., 
      seed=%eval(&seed.+&bs1_ebn.), simming=1);

    %blrm_ewoc(dose_list=&dose_list., ref_dose=&ref_dose., mcmcout=out_mcmc,
      dataout=next_dose, ud=&ud., od=&od., ewoc=&ewoc., debug=&debug.);

    /* EWOC records with batch label */
    data ewoc_temp; set next_dose; enroll_batch=&bs1_ebn.; run;
    data ewoc_record; set ewoc_record ewoc_temp; run;

    /* Dose selection (restrict candidates by skip rule, pick max TT_Mean among EWOC=TRUE) */
    proc sort data=dose_sort;  by dose; run;
    proc sort data=next_dose;  by dose; run;

    data tt1;
      merge dose_sort(in=a) next_dose(in=b);
      by dose;
      if a;
      if upcase(EWOC)='TRUE';
      %if &skip=0 %then %do; if ord <= &bs1_cur_ord.+1; %end;
      %else %do;              if ord <= &bs1_cur_ord.+2; %end;
    run;

    proc sort data=tt1; by descending TT_Mean; run;

    /* how many EWOC-eligible? */
    data _null_;
      if 0 then set tt1 nobs=nobs;
      call symputx('bs1_ewoc_left', nobs);
      stop;
    run;

    %if &bs1_ewoc_left = 0 %then %do;
      %let bs1_cur_dose = .;
      %let bs1_sample_expand = %eval(&bs1_sample_expand. + &cohortn.);
      %goto bs1_exit;
    %end;

    /* choose next */
    data tt2;
      set tt1;
      by descending TT_Mean;
      if _n_ = 1 then do;
        call symputx('bs1_cur_dose', dose);
        call symputx('bs1_cur_ord',  ord);
        output;
      end;
    run;

    /* counters */
    %let bs1_sample_expand = %eval(&bs1_sample_expand. + &cohortn.);
    %let bs1_ebn          = %eval(&bs1_ebn.+1);
    %let bs1_iter_count   = %eval(&bs1_iter_count.+1);

  %end; /* loop */

%bs1_exit:
  /* Final outputs */
  data dlt_record;  set dlt_record;  sim_batch=&sim_batch.;scenario="&scenario"; run;
  data ewoc_record; set ewoc_record; sim_batch=&sim_batch.;scenario="&scenario"; run;

  data sim_record;
  	scenario="&scenario";
    dose_select = &bs1_cur_dose.;
    sample_used = &bs1_sample_expand.;
    sim_batch   = &sim_batch.;
  run;

  /* cleanup temp tables (keep bad_p_events/dlt/ewoc/sim_record) */
  %if &debug. = 0 %then %do;
    proc datasets nolist nowarn;
      delete tox_gen tox_sum out_mcmc next_dose ewoc_temp tt1 tt2 dose_sort;
    quit;
  %end;

%mend blrm_sim_conti_one;


/*-----------------------------------------------*
 | Batch runner with timing summary and outputs  |
 *-----------------------------------------------*/
%macro blrm_sim_conti_n(
    /* How many simulation replicates */
    sim_total=100,

    /* ---- pass-through: Design & BLRM ---- */
    dose_list=%str(5,10,20,40,80),
    ref_dose=,
    mu1=, mu2=, v1=, v2=, rho=,
    alpha=0.05, nbi=1000, nmc=10000,
    ud=0.16, od=0.33, ewoc=0.25,

    /* ---- pass-through: Trial logistics ---- */
    start_dose=,
    max_sample=,
    cohortn=3,
    skip=0,
    seed=12345,

    /* ---- pass-through: truth via FCMP ---- */
    tox_fn=,
    scenario=Medium,

    /* ---- pass-through: jitter on actual dose ---- */
    disturb_actual_lo=-1, disturb_actual_hi=1,

    /* ---- pass-through: handling bad p ---- */
    abort_on_bad_p=0,

    /* ---- misc ---- */
    debug=0,
    cleanup=1,                 /* delete per-iteration temps from inner macro */

    /* ---- timing summary ---- */
    time_out=time_summary,     /* output table name for timing summary */
    time_sum=1                 /* 1=print total timing, 0=no print */
  );

  %local _i total_start total_end sim_start sim_end;

  /* Basic check */
  %if %sysevalf(%superq(sim_total)=,boolean) or %sysevalf(&sim_total <= 0) %then %do;
    %put ERROR: sim_total must be a positive integer.;
    %return;
  %end;

  /* Wipe old aggregate outputs if they exist */
  proc datasets nolist nowarn;
    delete dlt_record_all ewoc_record_all sim_record_all bad_p_events_all
           sim_summary sim_selection_dist &time_out time_record _time_row _time_total;
  quit;

  /* Initialize timing dataset with full schema (prevents APPEND warnings) */
  data time_record;
    length scenario $40;
    length sim_batch 8 start_dt 8 end_dt 8 elapsed_sec 8;
    format start_dt end_dt datetime19.;
    stop;
  run;

  /* Total timer start */
  %let total_start = %sysfunc(datetime());

  /*==== main simulation loop ====*/
  %do _i=1 %to &sim_total;

    %let sim_start = %sysfunc(datetime());

    %blrm_sim_conti_one(
      dose_list=&dose_list.,
      ref_dose=&ref_dose.,
      mu1=&mu1., mu2=&mu2., v1=&v1., v2=&v2., rho=&rho.,
      alpha=&alpha., nbi=&nbi., nmc=&nmc.,
      ud=&ud., od=&od., ewoc=&ewoc.,
      start_dose=&start_dose.,
      max_sample=&max_sample.,
      cohortn=&cohortn.,
      skip=&skip.,
      seed=&seed.,
      sim_batch=&_i.,

      tox_fn=&tox_fn.,
      scenario=&scenario.,

      disturb_actual_lo=&disturb_actual_lo.,
      disturb_actual_hi=&disturb_actual_hi.,

      abort_on_bad_p=&abort_on_bad_p.,
      debug=&debug.
    );

    %let sim_end = %sysfunc(datetime());

    /* record timing for this sim */
    data _time_row;
      length scenario $40;
      scenario = "&scenario";
      sim_batch = &_i;
      start_dt  = &sim_start;
      end_dt    = &sim_end;
      elapsed_sec = end_dt - start_dt;
      format start_dt end_dt datetime19.;
    run;
    proc append base=time_record data=_time_row force; run;

    /* initialize or append aggregates */
    %if &_i = 1 %then %do;
      data dlt_record_all;  set dlt_record;  run;
      data ewoc_record_all; set ewoc_record; run;
      data sim_record_all;  set sim_record;  run;

      %if %sysfunc(exist(bad_p_events)) %then %do;
        data bad_p_events_all; set bad_p_events; run;
      %end;
      %else %do;
        data bad_p_events_all; length scenario $40 msg $200; stop; run;
      %end;
    %end;
    %else %do;
      proc append base=dlt_record_all     data=dlt_record     force; run;
      proc append base=ewoc_record_all    data=ewoc_record    force; run;
      proc append base=sim_record_all     data=sim_record     force; run;
      %if %sysfunc(exist(bad_p_events)) %then %do;
        proc append base=bad_p_events_all data=bad_p_events   force; run;
      %end;
    %end;

    /* optional cleanup */
    %if &cleanup %then %do;
      proc datasets nolist nowarn;
        delete dose_sort tox_gen tox_sum out_mcmc next_dose ewoc_temp tt1 tt2 _time_row;
      quit;
    %end;
			%if &debug=0 %then %do;
			proc delete data=dlt_record ewoc_record sim_record;
			%end;
  %end; /* loop */

  /* Total timer end */
  %let total_end = %sysfunc(datetime());

  /* overall timing row */
  data _time_total;
    total_start = &total_start;
    total_end   = &total_end;
    total_elapsed_sec = total_end - total_start;
    avg_sec_per_sim   = total_elapsed_sec / &sim_total.;
    sims              = &sim_total.;
    format total_start total_end datetime19.;
  run;

  /* Combine overall + per-sim timing into &time_out */
  data &time_out;
    set _time_total time_record;
  run;

  /* Quick outcomes summary */
  proc sql;
    create table sim_summary as
    select 
      "&scenario" as scenario length=40,
      count(*)                                  as sims,
      sum(missing(dose_select))                 as n_early_stop,
      mean(sample_used)                         as avg_sample_used,
      calculated n_early_stop / calculated sims as early_stop_rate
    from sim_record_all;

    create table sim_selection_dist as
    select dose_select as dose,
           count(*) as n_select,
           calculated n_select /
             (select count(*) from sim_record_all where not missing(dose_select)) as pct_select
    from sim_record_all
    where not missing(dose_select)
    group by dose_select
    order by dose_select;
  quit;

  /* optional print: total timing only */
  %if &time_sum %then %do;
    title "Simulation Timing Summary (Total)";
    proc print data=&time_out noobs;
      where missing(sim_batch);  /* the single total row */
      var sims total_start total_end total_elapsed_sec avg_sec_per_sim;
    run;
    title;
  %end;

%mend blrm_sim_conti_n;

/*------------------*
 | Example (optional)
 *------------------*/
/*
%blrm_sim_conti_n(
  sim_total=50,
  dose_list=%str(5,10,20,40,80),
  ref_dose=20,
  mu1=-0.7, mu2=0, v1=%sysevalf(2**2), v2=1, rho=0,
  alpha=0.05, nbi=1000, nmc=10000,
  ud=0.18, od=0.30, ewoc=0.25,
  start_dose=10,
  max_sample=36,
  cohortn=3,
  skip=0,
  seed=12345,
  tox_fn=p_true,
  scenario=Medium,
  disturb_actual_lo=-1, disturb_actual_hi=1,
  abort_on_bad_p=1,
  debug=0,
  cleanup=1,
  time_out=time_summary,
  time_sum=1
);
*/
