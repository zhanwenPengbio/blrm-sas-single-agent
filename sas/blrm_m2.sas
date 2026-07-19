/*======================================================================*/
/* Core BLRM fitting macros with MCMC diagnostics                       */
/*                                                                      */
/* blrm_main                                                            */
/*   Single-chain PROC MCMC engine.                                     */
/*   diagnostics=1 saves within-chain diagnostics:                      */
/*     <diagout>_ess, <diagout>_mcse, <diagout>_autocorr,               */
/*     <diagout>_geweke.                                                */
/*                                                                      */
/* blrm_main_mc                                                         */
/*   Multi-chain wrapper for formal posterior inference.                */
/*   It calls blrm_main independently for each chain, pools posterior   */
/*   samples, saves chain-level diagnostics, and calculates classical   */
/*   Gelman-Rubin R-hat for log_alpha and log_beta.                     */
/*======================================================================*/

/*----------------------------------------------------------------------*/
/* blrm_main: single-chain BLRM posterior sampler                       */
/*                                                                      */
/* datain       : input data set with dose, n, and dltn                 */
/* dataout      : posterior-draw output data set                        */
/* ref_dose     : reference dose used to normalize dose                 */
/* mu1          : prior mean of log(alpha)                              */
/* mu2          : prior mean of log(beta)                               */
/* v1           : prior variance of log(alpha)                          */
/* v2           : prior variance of log(beta)                           */
/* rho          : prior covariance of log(alpha) and log(beta)          */
/* alpha        : tail probability for PROC MCMC posterior intervals    */
/* nbi          : burn-in iterations                                    */
/* nmc          : retained sampling iterations                          */
/* simming      : simulation indicator; suppresses standard ODS output  */
/* seed         : random-number seed                                    */
/* diagnostics  : 1=save within-chain diagnostics; 0=disable them       */
/* diagout      : diagnostic data-set prefix; if blank, use             */
/*                <dataout>_diag                                       */
/* diagplots    : 1=display trace and autocorrelation plots             */
/* ac_lags      : autocorrelation lags used in DIAG=AUTOCORR             */
/* init_random  : 1=initialize the chain randomly from the priors       */
/* postsumout   : posterior interval ODS data set when simming=0        */
/*----------------------------------------------------------------------*/
%macro blrm_main(
    datain=,
    dataout=,
    ref_dose=,
    mu1=,
    mu2=,
    v1=,
    v2=,
    rho=,
    alpha=0.05,
    nbi=2000,
    nmc=20000,
    simming=0,
    seed=9527,
    diagnostics=0,
    diagout=,
    diagplots=0,
    ac_lags=%str(1 5 10 50),
    init_random=0,
    postsumout=PostSumInt
);

    %local diag_ess diag_mcse diag_autocorr diag_geweke has_ods_output;

    %if %length(%superq(datain)) = 0 %then %do;
        %put ERROR: datain must be specified in blrm_main.;
        %return;
    %end;

    %if %length(%superq(dataout)) = 0 %then %do;
        %put ERROR: dataout must be specified in blrm_main.;
        %return;
    %end;

    %if %length(%superq(diagnostics)) = 0 %then %let diagnostics=0;
    %if %length(%superq(diagplots)) = 0 %then %let diagplots=0;
    %if %length(%superq(init_random)) = 0 %then %let init_random=0;

    %let has_ods_output=0;

    %if &diagnostics. = 1 %then %do;
        %if %length(%superq(diagout)) = 0 %then
            %let diagout=&dataout._diag;

        %let diag_ess      = &diagout._ess;
        %let diag_mcse     = &diagout._mcse;
        %let diag_autocorr = &diagout._autocorr;
        %let diag_geweke   = &diagout._geweke;
    %end;

    %if &diagplots. = 1 %then %do;
        ods graphics on;
    %end;

    /* Keep existing simulation behavior: suppress procedure output. */
    ods exclude all;

    %if &simming. = 0 %then %do;
        ods exclude close;
        %if %length(%superq(postsumout)) > 0 %then %do;
            ods output PostSumInt=&postsumout.;
            %let has_ods_output=1;
        %end;
    %end;

    %if &diagnostics. = 1 %then %do;
        %let has_ods_output=1;
        ods output
            ESS      = &diag_ess.
            MCSE     = &diag_mcse.
            AUTOCORR = &diag_autocorr.
            GEWEKE   = &diag_geweke.
        ;
    %end;

    proc mcmc
        data=&datain.
        stats(alpha=&alpha.)
        seed=&seed.
        nbi=&nbi.
        nmc=&nmc.
        outpost=&dataout.
        %if &init_random. = 1 %then init=random;
        %if &diagnostics. = 1 %then
            diag=(AUTOCORR(LAGS=(&ac_lags.)) ESS MCSE GEWEKE);
        %else
            diag=none;
        %if &diagplots. = 1 %then
            plots=(TRACE AUTOCORR);
        %else
            plots=none;
    ;

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

        delta = log_par[1] + exp(log_par[2]) * log(dose / &ref_dose.);
        pi_d = logistic(delta);

        model dltn ~ binomial(n,pi_d);
    run;

    %if &has_ods_output. = 1 %then %do;
        ods output close;
    %end;

    %if &diagplots. = 1 %then %do;
        ods graphics off;
    %end;

%mend blrm_main;


/*----------------------------------------------------------------------*/
/* blrm_main_mc: multi-chain BLRM analysis wrapper                      */
/*                                                                      */
/* datain         : input data set with dose, n, and dltn               */
/* dataout        : pooled posterior-draw output data set               */
/* ref_dose       : reference dose used to normalize dose                */
/* mu1, mu2       : prior means of log(alpha), log(beta)                 */
/* v1, v2         : prior variances of log(alpha), log(beta)             */
/* rho            : prior covariance of log(alpha) and log(beta)         */
/* alpha          : tail probability for posterior intervals             */
/* nbi            : burn-in iterations per chain                         */
/* nmc            : retained sampling iterations per chain               */
/* nchain         : number of independent PROC MCMC runs                 */
/* simming        : simulation indicator; normally use 0 for this macro  */
/* seed           : seed for chain 1; subsequent chains use seed+1, ...  */
/* diagnostics    : 1=save single-chain diagnostics for every chain      */
/* chain_diagout  : prefix for pooled chain-level diagnostic data sets    */
/* diagplots      : reserved for single-chain inspection; keep 0 here     */
/* ac_lags        : autocorrelation lags used in each chain              */
/* init_random    : 1=randomize starting values independently by chain   */
/* rhat_cutoff    : R-hat threshold for the convergence flag             */
/* postsumout     : pooled posterior parameter summary                    */
/* diagout        : R-hat output data set                                 */
/* debug          : 1=retain internal per-chain temporary data sets       */
/*                                                                      */
/* Additional outputs when diagnostics=1:                               */
/*   <chain_diagout>_ess                                                 */
/*   <chain_diagout>_mcse                                                */
/*   <chain_diagout>_autocorr                                            */
/*   <chain_diagout>_geweke                                              */
/*----------------------------------------------------------------------*/
%macro blrm_main_mc(
    datain=,
    dataout=,
    ref_dose=,
    mu1=,
    mu2=,
    v1=,
    v2=,
    rho=,
    alpha=0.05,
    nbi=2000,
    nmc=20000,
    nchain=4,
    simming=0,
    seed=9527,
    diagnostics=1,
    chain_diagout=,
    diagplots=0,
    ac_lags=%str(1 5 10 50),
    init_random=1,
    rhat_cutoff=1.01,
    postsumout=,
    diagout=,
    debug=0
);

    %local i chain_seed chain_list internal_simming
           chain_ess_list chain_mcse_list chain_autocorr_list chain_geweke_list;

    %if %length(%superq(postsumout)) = 0 %then
        %let postsumout=&dataout._postsum;

    %if %length(%superq(diagout)) = 0 %then
        %let diagout=&dataout._rhat;

    %if %length(%superq(chain_diagout)) = 0 %then
        %let chain_diagout=&dataout._chain_diag;

    %if %length(%superq(datain)) = 0 %then %do;
        %put ERROR: datain must be specified in blrm_main_mc.;
        %return;
    %end;

    %if %length(%superq(dataout)) = 0 %then %do;
        %put ERROR: dataout must be specified in blrm_main_mc.;
        %return;
    %end;

    %if %length(%superq(nchain)) = 0 %then %let nchain=1;

    %if %sysevalf(&nchain. < 1) %then %do;
        %put ERROR: nchain must be greater than or equal to 1.;
        %return;
    %end;

    %let chain_list=;
    %let chain_ess_list=;
    %let chain_mcse_list=;
    %let chain_autocorr_list=;
    %let chain_geweke_list=;

    /* Show plots only when explicitly requested. */
    %if &diagplots. = 1 %then %let internal_simming=0;
    %else %let internal_simming=1;

    /*---------------------------------------------------------------*/
    /* Run independent PROC MCMC chains.                             */
    /* simming=1 is used internally to avoid one ODS report per run. */
    /*---------------------------------------------------------------*/
    %do i=1 %to &nchain.;

        %let chain_seed=%eval(&seed. + &i. - 1);

        %put NOTE: ======================================================;
        %put NOTE: Running BLRM chain &i. of &nchain. with seed=&chain_seed.;
        %put NOTE: ======================================================;

        %blrm_main(
            datain=&datain.,
            dataout=work._blrm_mc_chain&i.,
            ref_dose=&ref_dose.,
            mu1=&mu1.,
            mu2=&mu2.,
            v1=&v1.,
            v2=&v2.,
            rho=&rho.,
            alpha=&alpha.,
            nbi=&nbi.,
            nmc=&nmc.,
            simming=&internal_simming.,
            seed=&chain_seed.,
            diagnostics=&diagnostics.,
            diagout=work._blrm_mc_diag&i.,
            diagplots=&diagplots.,
            ac_lags=&ac_lags.,
            init_random=&init_random.,
            postsumout=
        );

        data work._blrm_mc_chain&i.;
            set work._blrm_mc_chain&i.;
            chain=&i.;
            chain_seed=&chain_seed.;
            chain_iteration=_n_;
        run;

        %let chain_list=&chain_list. work._blrm_mc_chain&i.;

        %if &diagnostics. = 1 %then %do;

            data work._blrm_mc_ess&i.;
                set work._blrm_mc_diag&i._ess;
                chain=&i.;
                chain_seed=&chain_seed.;
            run;

            data work._blrm_mc_mcse&i.;
                set work._blrm_mc_diag&i._mcse;
                chain=&i.;
                chain_seed=&chain_seed.;
            run;

            data work._blrm_mc_autocorr&i.;
                set work._blrm_mc_diag&i._autocorr;
                chain=&i.;
                chain_seed=&chain_seed.;
            run;

            data work._blrm_mc_geweke&i.;
                set work._blrm_mc_diag&i._geweke;
                chain=&i.;
                chain_seed=&chain_seed.;
            run;

            %let chain_ess_list=&chain_ess_list. work._blrm_mc_ess&i.;
            %let chain_mcse_list=&chain_mcse_list. work._blrm_mc_mcse&i.;
            %let chain_autocorr_list=&chain_autocorr_list. work._blrm_mc_autocorr&i.;
            %let chain_geweke_list=&chain_geweke_list. work._blrm_mc_geweke&i.;

        %end;

    %end;

    /*---------------------------------------------------------------*/
    /* Pool posterior draws from all chains.                         */
    /*---------------------------------------------------------------*/
    data &dataout.;
        set &chain_list.;
    run;

    proc sort data=&dataout.;
        by chain chain_iteration;
    run;

    /*---------------------------------------------------------------*/
    /* Pool raw within-chain diagnostics, retaining chain labels.    */
    /*---------------------------------------------------------------*/
    %if &diagnostics. = 1 %then %do;

        data &chain_diagout._ess;
            set &chain_ess_list.;
        run;

        data &chain_diagout._mcse;
            set &chain_mcse_list.;
        run;

        data &chain_diagout._autocorr;
            set &chain_autocorr_list.;
        run;

        data &chain_diagout._geweke;
            set &chain_geweke_list.;
        run;

    %end;

    /*---------------------------------------------------------------*/
    /* Convert posterior samples to long form.                       */
    /*---------------------------------------------------------------*/
    data work._blrm_mc_long;
        set &dataout.;

        length parameter $12;

        parameter='log_alpha';
        value=log_alpha;
        output;

        parameter='log_beta';
        value=log_beta;
        output;

        keep chain chain_seed chain_iteration parameter value;
    run;

    /*---------------------------------------------------------------*/
    /* Pooled posterior parameter summaries.                         */
    /*---------------------------------------------------------------*/
    proc sort data=work._blrm_mc_long;
        by parameter;
    run;

    proc means data=work._blrm_mc_long noprint;
        by parameter;
        var value;
        output out=&postsumout.(drop=_TYPE_ _FREQ_)
            n=N_Draws
            mean=Mean
            std=StdDev
            min=Min
            q1=Q1
            median=Median
            q3=Q3
            max=Max;
    run;

    data &postsumout.;
        set &postsumout.;

        label
            parameter='Posterior Parameter'
            N_Draws='Number of Pooled Draws'
            Mean='Posterior Mean'
            StdDev='Posterior SD'
            Min='Posterior Minimum'
            Q1='Posterior Q1'
            Median='Posterior Median'
            Q3='Posterior Q3'
            Max='Posterior Maximum';
    run;

    /*---------------------------------------------------------------*/
    /* Classical Gelman-Rubin R-hat.                                */
    /* W       = mean within-chain variance                         */
    /* B       = n * variance of chain means                        */
    /* VarPlus = ((n-1)/n) * W + B/n                                */
    /* Rhat    = max(1, sqrt(VarPlus/W))                            */
    /*---------------------------------------------------------------*/
    proc sort data=work._blrm_mc_long;
        by parameter chain;
    run;

    proc means data=work._blrm_mc_long noprint;
        by parameter chain;
        var value;
        output out=work._blrm_mc_chain_stats(drop=_TYPE_ _FREQ_)
            n=n_draw
            mean=chain_mean
            var=chain_var;
    run;

    proc sort data=work._blrm_mc_chain_stats;
        by parameter;
    run;

    proc means data=work._blrm_mc_chain_stats noprint;
        by parameter;
        var n_draw chain_var chain_mean;
        output out=work._blrm_mc_rhat_base(drop=_TYPE_ _FREQ_)
            n(chain_mean)=N_Chain
            min(n_draw)=N_Draw_Min
            max(n_draw)=N_Draw_Max
            mean(chain_var)=W
            var(chain_mean)=Var_Chain_Mean;
    run;

    data &diagout.;
        set work._blrm_mc_rhat_base;

        length Rhat_Status $48 Converged $3;

        Rhat_Cutoff=&rhat_cutoff.;

        if N_Chain < 2 then do;
            N_Draw_Per_Chain=N_Draw_Min;
            Between_Chain_B=.;
            Var_Plus=.;
            Rhat=.;
            Converged='NA';
            Rhat_Status='Not available: nchain < 2';
        end;
        else if N_Draw_Min ne N_Draw_Max then do;
            N_Draw_Per_Chain=.;
            Between_Chain_B=.;
            Var_Plus=.;
            Rhat=.;
            Converged='No';
            Rhat_Status='Unavailable: unequal chain lengths';
        end;
        else if missing(W) or W <= 0 then do;
            N_Draw_Per_Chain=N_Draw_Min;
            Between_Chain_B=.;
            Var_Plus=.;
            Rhat=.;
            Converged='No';
            Rhat_Status='Unavailable: invalid within-chain variance';
        end;
        else do;
            N_Draw_Per_Chain=N_Draw_Min;
            Between_Chain_B=N_Draw_Per_Chain * Var_Chain_Mean;
            Var_Plus=((N_Draw_Per_Chain - 1) / N_Draw_Per_Chain) * W
                + Between_Chain_B / N_Draw_Per_Chain;
            Rhat=max(1, sqrt(Var_Plus / W));

            if Rhat <= Rhat_Cutoff then do;
                Converged='Yes';
                Rhat_Status='Rhat within cutoff';
            end;
            else do;
                Converged='No';
                Rhat_Status='Rhat exceeds cutoff';
            end;
        end;

        label
            parameter='Posterior Parameter'
            N_Chain='Number of Chains'
            N_Draw_Min='Minimum Draws per Chain'
            N_Draw_Max='Maximum Draws per Chain'
            N_Draw_Per_Chain='Draws per Chain Used in Rhat'
            W='Within-Chain Variance'
            Var_Chain_Mean='Variance of Chain Means'
            Between_Chain_B='Between-Chain Variance Component'
            Var_Plus='Estimated Marginal Posterior Variance'
            Rhat='Gelman-Rubin Rhat'
            Rhat_Cutoff='Rhat Convergence Cutoff'
            Converged='Rhat Convergence Flag'
            Rhat_Status='Rhat Diagnostic Status';

        keep
            parameter
            N_Chain
            N_Draw_Min
            N_Draw_Max
            N_Draw_Per_Chain
            W
            Var_Chain_Mean
            Between_Chain_B
            Var_Plus
            Rhat
            Rhat_Cutoff
            Converged
            Rhat_Status;
    run;

    /* Restore visible ODS output after formal analysis. */
    %if &simming. = 0 %then %do;
        ods exclude none;
    %end;

    /* Retain only user-facing outputs unless debug=1. */
    %if &debug. = 0 %then %do;
        proc datasets lib=work nolist;
            delete
                %do i=1 %to &nchain.;
                    _blrm_mc_chain&i.
                    _blrm_mc_diag&i._ess
                    _blrm_mc_diag&i._mcse
                    _blrm_mc_diag&i._autocorr
                    _blrm_mc_diag&i._geweke
                    _blrm_mc_ess&i.
                    _blrm_mc_mcse&i.
                    _blrm_mc_autocorr&i.
                    _blrm_mc_geweke&i.
                %end;
                _blrm_mc_long
                _blrm_mc_chain_stats
                _blrm_mc_rhat_base;
        quit;
    %end;

%mend blrm_main_mc;


%macro loop(values);    
                                                                                                                
     /* Count the number of values in the string */                                                                                                                                   
     %let count=%sysfunc(countw(&values)); 

     /* Loop through the total number of values */                                                                                         
     %do i = 1 %to &count;                                                                                                              
      %let value=%qscan(&values,&i,%str(,));                                                                                            
      %put &value;                                                                                                                      
     %end;                                                                                                                              
                                                                                                                                        
%mend;  


%macro make_dose_map(dose_list=, out=work._dose_map);
    %local n_dose i;
    %let n_dose = %sysfunc(countw(%superq(dose_list), %str(,)));

    data &out.;
        length dose_raw $32;
        %do i = 1 %to &n_dose.;
            dose_id = &i.;
            dose_raw = strip("%qscan(%superq(dose_list), &i., %str(,))");
            dose = input(dose_raw, ?? best32.);
            if missing(dose) then do;
                put "ERROR: Invalid dose value in dose_list: " dose_raw=;
                stop;
            end;
            output;
        %end;
        drop dose_raw;
    run;
%mend;


%macro make_dose_truth_map(dose_list=, true_tox=, out=work._dose_truth_map);
    %make_dose_map(dose_list=&dose_list., out=work._dose_map_base);

    proc sort data=work._dose_map_base;
        by dose;
    run;

    proc sort data=&true_tox.(keep=dose true_tox)
              out=work._true_tox_srt nodupkey;
        by dose;
    run;

    data &out.;
        merge work._dose_map_base(in=a)
              work._true_tox_srt(in=b);
        by dose;
        if a;
    run;
%mend;





/*blrm_stat: a macro that summarize posterior dlt from blrm_main.*/
/*dose_list: the dose levels that is wanted to summarize the posterior dlt, not necessary the same dose levels in blrm_main.*/
/*(example of dose_list: %str(2,4,5,12))*/
/*mcmcout: the output dataset from*/
/*dataout: the summary of blrm output*/
/*ref_dose: reference dose, for scaling the dose range.*/
/*stat: the statistics want to summarize*/
/*debug: if set to 0, temporary dataset will be deleted.*/


%macro blrm_stat(
    dose_list=%str(),
    ref_dose=,
    mcmcout=,
    stat=%str(mean, std, q1, median, q3),
    dataout=,
    debug=0
);


    %make_dose_map(dose_list=&dose_list., out=work._dose_map_stat);

    proc sql;
        create table work._temp_long_stat as
        select
            a.Iteration,
            b.dose_id,
            b.dose,
            logistic(a.log_alpha + exp(a.log_beta) * log(b.dose / &ref_dose.)) as pi_d
        from &mcmcout. as a,
             work._dose_map_stat as b
        ;
    quit;

    proc sort data=work._temp_long_stat;
        by dose_id dose;
    run;

    proc means data=work._temp_long_stat noprint;
        by dose_id dose;
        var pi_d;
        output out=&dataout.(drop=_TYPE_ _FREQ_)
            mean   = Mean
            std    = StdDev
            q1     = Q1
            median = Median
            q3     = Q3
        ;
    run;

    %if &debug.=0 %then %do;
        proc datasets lib=work nolist;
            delete _dose_map_stat _temp_long_stat;
        quit;
    %end;

%mend;

/*blrm_ewoc: a macro that summarize under-dose, target dose, over dose.*/
/*dose_list: the dose levels that is wanted to summarize the posterior dlt, not necessary the same dose levels in blrm_main.*/
/*(example of dose_list: %str(2,4,5,12))*/
/*mcmcout: the output dataset from*/
/*dataout: the summary of blrm output*/
/*ud: under dosing boundary*/
/*od: over dosing boundary*/
/*ewoc: the escalation with overdose control RULE, control the posterior dlt not exceed some specific value*/
/*debug: if set to 0, temporary dataset will be deleted.*/

/*blrm_ewoc: a macro that summarize under-dose, target dose, over dose.*/
/*blrm_ewoc: a macro that summarize under-dose, target dose, over dose.*/
%macro blrm_ewoc(
    dose_list=,
    ref_dose=,
    mcmcout=out_mcmc,
    dataout=d1,
    ud=0.16,
    od=0.33,
    ewoc=0.25,
    debug=0
);

    %make_dose_map(dose_list=&dose_list., out=work._dose_map_ewoc);

    proc sql;
        create table work._ewoc_long as
        select
            a.Iteration,
            b.dose_id,
            b.dose,
            logistic(a.log_alpha + exp(a.log_beta) * log(b.dose / &ref_dose.)) as pi_d
        from &mcmcout. as a,
             work._dose_map_ewoc as b
        ;
    quit;

    data work._ewoc_long2;
        set work._ewoc_long;
        ud_ind = (pi_d < &ud.);
        tt_ind = (pi_d >= &ud. and pi_d < &od.);
        od_ind = (pi_d >= &od.);
    run;

    proc sort data=work._ewoc_long2;
        by dose_id dose;
    run;

    proc means data=work._ewoc_long2 noprint;
        by dose_id dose;
        var ud_ind tt_ind od_ind;
        output out=work._ewoc_sum
            mean(ud_ind)=UD_Mean
            mean(tt_ind)=TT_Mean
            mean(od_ind)=OD_Mean
        ;
    run;

    data &dataout.;
        set work._ewoc_sum;
        length EWOC $5;
        if OD_Mean <= &ewoc. then EWOC = 'True';
        else EWOC = 'False';
		check1 = sum(UD_Mean, TT_Mean, OD_Mean);
        keep dose_id  dose UD_Mean TT_Mean OD_Mean EWOC %if &debug.=1 %then check1; ;
        %if &debug.=1 %then %do;
            check1 = sum(UD_Mean, TT_Mean, OD_Mean);
        %end;
    run;

    proc sort data=&dataout.;
        by dose;
    run;

    %if &debug.=0 %then %do;
        proc datasets lib=work nolist;
            delete _dose_map_ewoc _ewoc_long _ewoc_long2 _ewoc_sum;
        quit;
    %end;

%mend;


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

    %local bs1_ebn bs1_iter_count bs1_max_iter bs1_sample_expand current_cohortn
           bs1_cur_dose bs1_cur_ord bs1_cur_tox bs1_ewoc_left
           bs1_next_ord bs1_next_dose _nmiss;

    %let bs1_ebn = 1;
    %let bs1_sample_expand = 0;
    %let bs1_max_iter = %eval(&max_sample./&cohortn. + 20);
    %let bs1_iter_count = 1;
    %let current_cohortn = &cohortn.;
	data dlt_record;
    	stop;
	run;

	data ewoc_record;
	    stop;
	run;

	data sim_record;
	    stop;
	run;
    %make_dose_truth_map(
        dose_list=&dose_list.,
        true_tox=&true_tox.,
        out=work._dose_truth_map
    );

    proc sql noprint;
        select count(*) into :_nmiss trimmed
        from work._dose_truth_map
        where missing(true_tox);
    quit;

    %if %eval(&_nmiss. > 0) %then %do;
        %put ERROR: true_tox is missing for one or more doses in dose_list.;
        proc print data=work._dose_truth_map;
            where missing(true_tox);
        run;
        %return;
    %end;

    proc sql noprint;
        select dose_id, dose, true_tox
        into :bs1_cur_ord trimmed,
             :bs1_cur_dose trimmed,
             :bs1_cur_tox trimmed
        from work._dose_truth_map
        where abs(dose - &start_dose.) < 1e-12;
    quit;

    %if %superq(bs1_cur_ord)= %then %do;
        %put ERROR: start_dose=&start_dose. is not found in dose_list/true_tox.;
        %return;
    %end;

    data dlt_record;
        stop;
    run;

    data ewoc_record;
        stop;
    run;

    %do %until (&bs1_sample_expand.>=&max_sample. or &bs1_iter_count. > &bs1_max_iter.);

        %let current_cohortn = &cohortn.;
        %if %eval(&max_sample.-&bs1_sample_expand.)<=&cohortn. and %eval(&max_sample.-&bs1_sample_expand.)>0 %then %do;
            %let current_cohortn = %eval(&max_sample.-&bs1_sample_expand.);
        %end;

        data tox_gen;
            dose = &bs1_cur_dose.;
            n = &current_cohortn.;
            call streaminit(%eval(&seed. + &sim_batch.*10000 + &bs1_ebn.*100));
            dltn = rand("BINOMIAL", &bs1_cur_tox., n);
            true_tox = &bs1_cur_tox.;
            enroll_batch = &bs1_ebn.;
        run;

        data dlt_record;
            set dlt_record tox_gen;
        run;

        proc sql;
            create table tox_sum as
            select dose, sum(n) as n, sum(dltn) as dltn
            from dlt_record
            group by dose;
        quit;

        %blrm_main(
            datain=tox_sum,
            dataout=out_mcmc,
            ref_dose=&ref_dose.,
            mu1=&mu1., mu2=&mu2., v1=&v1., v2=&v2., rho=&rho.,
            alpha=&alpha., nbi=&nbi., nmc=&nmc.,
            seed=%eval(&seed.+&bs1_ebn.), simming=1
        );

        %blrm_ewoc(
            dose_list=&dose_list.,
            ref_dose=&ref_dose.,
            mcmcout=out_mcmc,
            dataout=next_dose,
            ud=&ud.,
            od=&od.,
            ewoc=&ewoc.,
            debug=&debug.
        );

        data ewoc_temp;
            set next_dose;
            enroll_batch = &bs1_ebn.;
        run;

        data ewoc_record;
            set ewoc_record ewoc_temp;
        run;

        data work._ewoc_pick;
            set next_dose;
            where EWOC = 'True';
            %if &skip.=0 %then %do;
                if dose_id <= (&bs1_cur_ord. + 1);
            %end;
            %else %if &skip.=1 %then %do;
                if dose_id <= (&bs1_cur_ord. + 2);
            %end;
        run;

        proc sort data=work._ewoc_pick;
            by descending TT_Mean descending dose_id;
        run;

        data _null_;
            if 0 then set work._ewoc_pick nobs=nobs;
            call symputx('bs1_ewoc_left', nobs);
            stop;
        run;

        %if &bs1_ewoc_left. = 0 %then %do;
            %let bs1_cur_dose = .;
            %let bs1_sample_expand = %eval(&bs1_sample_expand. + &current_cohortn.);
            %goto bs1_exit;
        %end;

        data _null_;
            set work._ewoc_pick(obs=1);
            call symputx('bs1_next_ord', dose_id);
            call symputx('bs1_next_dose', dose);
        run;

        proc sql noprint;
            select true_tox
            into :bs1_cur_tox trimmed
            from work._dose_truth_map
            where dose_id = &bs1_next_ord.;
        quit;

        %let bs1_cur_ord  = &bs1_next_ord.;
        %let bs1_cur_dose = &bs1_next_dose.;

        %let bs1_sample_expand = %eval(&bs1_sample_expand. + &current_cohortn.);
        %let bs1_ebn = %eval(&bs1_ebn.+1);
        %let bs1_iter_count = %eval(&bs1_iter_count.+1);

    %end;

%bs1_exit:

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
        proc datasets lib=work nolist;
            delete _dose_map_base _true_tox_srt _dose_truth_map _ewoc_pick
                   tox_gen tox_sum out_mcmc next_dose ewoc_temp;
        quit;
    %end;

%mend blrm_sim_one;


%macro blrm_sim_n(dose_list=%str(),ref_dose=,true_tox=,start_dose=,max_sample=,cohortn=,skip=,
        mu1=, mu2=, v1=, v2=, rho= ,
        alpha=0.05, nbi=1000, nmc = 10000, 
        ud=0.16,od=0.33,ewoc=0.25,
        seed=,sim_time=4,debug=0, time_out=time_summary,time_sum=1);
        
    %local bs_n_i start_time end_time total_time;
    data dlt_all;
        stop;
    run;
    data ewoc_all;
        stop;
    run;
    data sim_all;
        stop;
    run;
    data &time_out;
        stop;
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
    dlt_data=,              /* Cohort-level data from blrm_sim_n (dlt_all) */
    true_tox=,              /* Dataset with true toxicity probabilities */
    target_range_low=0.16,  /* Lower bound of target toxicity range */
    target_range_high=0.33, /* Upper bound of target toxicity range; also overdose threshold */
    acceptable_low=0.05,    /* Lower bound for an acceptable, non-over-toxic dose */
    out_summary=,           /* Output dataset name */
    debug=0
);

    %local _n_dup_truth _n_bad_truth _n_unmapped_select _n_unmapped_exposure;

    /*
      Category convention is deliberately identical to blrm_ewoc:
        under-dose:          pi < target_range_low
        target dose:         target_range_low <= pi < target_range_high
        acceptable dose:     acceptable_low   <= pi < target_range_high
        over-toxic dose:     pi >= target_range_high
    */

    /* Step 0: Basic truth-data checks.  A duplicate dose would duplicate rows in SQL joins. */
    proc sql noprint;
        select count(*) - count(distinct dose)
          into :_n_dup_truth trimmed
        from &true_tox.;

        select count(*)
          into :_n_bad_truth trimmed
        from &true_tox.
        where missing(dose)
           or missing(true_tox)
           or true_tox < 0
           or true_tox > 1;
    quit;

    %if &_n_dup_truth. > 0 %then %do;
        %put ERROR: blrm_sim_sum: true_tox contains duplicate dose values.;
        %return;
    %end;

    %if &_n_bad_truth. > 0 %then %do;
        %put ERROR: blrm_sim_sum: true_tox must contain nonmissing dose values and probabilities in [0,1].;
        %return;
    %end;

    /* Step 1: Define dose categories using the same half-open intervals as EWOC. */
    data work._dose_categories;
        set &true_tox.(keep=dose true_tox);

        is_target      = (true_tox >= &target_range_low. and true_tox < &target_range_high.);
        is_acceptable  = (true_tox >= &acceptable_low.   and true_tox < &target_range_high.);
        is_over_toxic  = (true_tox >= &target_range_high.);
    run;

    /* One final recommendation record per simulated trial is expected from blrm_sim_n. */
    proc sql;
        create table work._sim_final as
        select sim_batch,
               dose_select,
               sample_used
        from &sim_data.;
    quit;

    /* Step 2: Classify the final recommendation.
       A missing dose_select means that no dose was recommended; it must contribute 0,
       rather than missing, to the selection probabilities so every simulated trial stays
       in the denominator. */
    proc sql;
        create table work._dose_selection_raw as
        select a.sim_batch,
               a.dose_select,
               a.sample_used,
               b.is_target,
               b.is_acceptable,
               b.is_over_toxic
        from work._sim_final as a
        left join work._dose_categories as b
          on a.dose_select = b.dose;
    quit;

    proc sql noprint;
        select count(*)
          into :_n_unmapped_select trimmed
        from work._dose_selection_raw
        where not missing(dose_select)
          and missing(is_target);
    quit;

    %if &_n_unmapped_select. > 0 %then %do;
        %put ERROR: blrm_sim_sum: one or more nonmissing dose_select values are absent from true_tox.;
        %return;
    %end;

    data work._dose_selection;
        set work._dose_selection_raw;

        no_dose_selected = missing(dose_select);

        if no_dose_selected then do;
            is_target     = 0;
            is_acceptable = 0;
            is_over_toxic = 0;
        end;
    run;

    proc means data=work._dose_selection noprint;
        var is_target is_acceptable is_over_toxic no_dose_selected sample_used;
        output out=work._selection_summary
            n(is_target)               = N_Simulations
            mean(is_target)            = Prop_Correct_Dose
            mean(is_acceptable)        = Prop_Acceptable_Dose
            mean(is_over_toxic)        = Prop_Over_Toxic_Dose
            mean(no_dose_selected)     = Prop_No_Dose_Selected
            mean(sample_used)          = Avg_Sample_Used;
    run;

    /* Step 3: Calculate patient exposure to truly over-toxic doses.
       Left-joining back to all simulated trials assigns 0 exposure to any trial without
       a record in dlt_data, preserving the all-trial denominator. */
    proc sql noprint;
        select count(*)
          into :_n_unmapped_exposure trimmed
        from &dlt_data. as d
        left join work._dose_categories as c
          on d.dose = c.dose
        where missing(c.dose);
    quit;

    %if &_n_unmapped_exposure. > 0 %then %do;
        %put ERROR: blrm_sim_sum: one or more doses in dlt_data are absent from true_tox.;
        %return;
    %end;

    proc sql;
        create table work._over_toxic_exposure_raw as
        select d.sim_batch,
               sum(case when c.is_over_toxic = 1 then d.n else 0 end) as n_over_toxic
        from &dlt_data. as d
        left join work._dose_categories as c
          on d.dose = c.dose
        group by d.sim_batch;

        create table work._over_toxic_exposure as
        select s.sim_batch,
               coalesce(e.n_over_toxic, 0) as n_over_toxic
        from work._sim_final as s
        left join work._over_toxic_exposure_raw as e
          on s.sim_batch = e.sim_batch;
    quit;

    proc means data=work._over_toxic_exposure noprint;
        var n_over_toxic;
        output out=work._exposure_summary
            mean = Avg_Patients_Over_Toxic;
    run;

    /* Step 4: Return a long-format OC summary. */
    data work._selection_metrics;
        length Metric $50 Value 8;
        set work._selection_summary(keep=N_Simulations
                                         Prop_Correct_Dose
                                         Prop_Acceptable_Dose
                                         Prop_Over_Toxic_Dose
                                         Prop_No_Dose_Selected
                                         Avg_Sample_Used);

        Metric = 'N_Simulations';            Value = N_Simulations;            output;
        Metric = 'Prop_Correct_Dose';        Value = Prop_Correct_Dose;        output;
        Metric = 'Prop_Acceptable_Dose';     Value = Prop_Acceptable_Dose;     output;
        Metric = 'Prop_Over_Toxic_Dose';     Value = Prop_Over_Toxic_Dose;     output;
        Metric = 'Prop_No_Dose_Selected';    Value = Prop_No_Dose_Selected;    output;
        Metric = 'Avg_Sample_Used';          Value = Avg_Sample_Used;          output;

        keep Metric Value;
    run;

    data work._exposure_metric;
        length Metric $50 Value 8;
        set work._exposure_summary(keep=Avg_Patients_Over_Toxic);
        Metric = 'Avg_Patients_Over_Toxic';
        Value  = Avg_Patients_Over_Toxic;
        keep Metric Value;
    run;

    data &out_summary.;
        length Parameters $200;
        set work._selection_metrics work._exposure_metric;

        Parameters = cats(
            'Acceptable: [', "&acceptable_low.", ', ', "&target_range_high.",
            ') | Target: [', "&target_range_low.", ', ', "&target_range_high.",
            ') | Over: [', "&target_range_high.", ', 1]'
        );
    run;

    /* Step 5: Clean up. */
    %if &debug. = 0 %then %do;
        proc datasets library=work nolist;
            delete _dose_categories _sim_final _dose_selection_raw _dose_selection
                   _selection_summary _over_toxic_exposure_raw _over_toxic_exposure
                   _exposure_summary _selection_metrics _exposure_metric;
        quit;
    %end;

%mend blrm_sim_sum;

