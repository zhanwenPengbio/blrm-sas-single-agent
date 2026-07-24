options obs=100;

/* Sample BLRM pipeline outputs, shaped like the datasets the export
   utility in sas/blrm_export.sas is designed to write out. Dose / n / dltn
   follow the input structure documented in the repository README. These are
   small toy values so the bundle is self-contained. */

/* INFO: the operating-characteristics summary (blrm_sim_sum output shape) */
data info;
  length Metric $40 Value 8 Parameters $60;
  Parameters = "Acceptable: >=0.05 | Target: 0.16-0.33 | Over: >0.33";
  Metric = "Prop_Correct_Dose";       Value = 0.62; output;
  Metric = "Prop_Acceptable_Dose";    Value = 0.88; output;
  Metric = "Prop_Over_Toxic_Dose";    Value = 0.07; output;
  Metric = "Avg_Patients_Over_Toxic"; Value = 1.30; output;
run;

/* DLT: patient-level enrollment record across simulated trials */
data dlt_all;
  input sim_batch dose n dltn;
  datalines;
1 5  3 0
1 10 6 1
2 5  3 0
2 10 6 1
2 20 6 2
3 5  3 0
3 10 6 1
;
run;

/* EWOC: under / target / over-dose summary per dose */
data ewoc_all;
  input sim_batch dose ud_mean tt_mean od_mean;
  length EWOC $5;
  EWOC = "True";
  datalines;
1 5  0.71 0.27 0.02
1 10 0.34 0.55 0.11
2 5  0.68 0.29 0.03
2 20 0.10 0.48 0.42
;
run;

/* SIM: one row per simulated trial with the selected dose */
data sim_all;
  input sim_batch dose_select sample_used;
  datalines;
1 10 24
2 20 27
3 10 21
;
run;

/* DOSE_SELECTION: how often each dose was chosen */
data dose_selection;
  input dose n_select pct_select;
  datalines;
10 2 0.667
20 1 0.333
;
run;

/* SIM_RESULTS: headline metrics for the run */
data sim_results;
  length scenario $20;
  input scenario $ sims avg_sample_used early_stop_rate;
  datalines;
Medium 3 24.0 0.0
;
run;
