# BLRM SAS Macro User Guide

This guide describes the current single-agent Bayesian logistic regression
model (BLRM) workflow in `sas/blrm_m2.sas` and the Excel export utility in
`sas/blrm_export.sas`.

The implementation supports two distinct activities:

- formal interim analysis using multiple independent MCMC chains and
  convergence diagnostics; and
- simulation-based design evaluation using a faster single-chain engine within
  each virtual trial.

The numerical examples are illustrative. They are not design recommendations.

## 1. Software requirements

- SAS 9.4 or later
- `PROC MCMC`
- ODS diagnostic tables when ESS, MCSE, autocorrelation, or Geweke diagnostics
  are requested
- the XLSX engine when `result_export` is used

Load the macro library at the beginning of a SAS session:

```sas
%include "path/to/sas/blrm_m2.sas";
%include "path/to/sas/blrm_export.sas";
```

## 2. Model and parameterization

For dose `d` and reference dose `d_ref`, the fitted model is

```text
logit(pi_d) = log_alpha + exp(log_beta) * log(d / d_ref)
```

with prior

```text
(log_alpha, log_beta)' ~ MVN(mu, Sigma)

Sigma = [ v1   rho ]
        [ rho  v2  ]
```

`v1` and `v2` are prior variances. `rho` is the prior covariance, not a
correlation coefficient. The covariance matrix must be positive definite.

The reference dose is a scaling choice and does not need to be an observed
dose. It must be positive, as must every modeled dose.

## 3. Input data

The fitting macros use one row per evaluated dose level.

| Variable | Type | Description |
|---|---|---|
| `dose` | Numeric | Dose level |
| `n` | Numeric | Number of treated participants |
| `dltn` | Numeric | Number of participants with a DLT |

```sas
data observed;
    input dose n dltn;
    datalines;
2 6 1
4 3 0
;
run;
```

Use cumulative DLT data at each interim analysis. The model is refitted from
the original prior using all currently available observations.

## 4. Macro roles

| Macro | Role | Recommended use |
|---|---|---|
| `blrm_main` | One `PROC MCMC` chain | Simple analysis, simulation internals, or detailed inspection of one chain |
| `blrm_main_mc` | Independent chains, pooled draws, posterior summary, classical Gelman-Rubin R-hat | Formal interim analysis |
| `blrm_stat` | Posterior toxicity summaries at candidate doses | Descriptive posterior review |
| `blrm_ewoc` | Underdosing, target-toxicity, and overdosing probabilities plus EWOC classification | Dose-decision support |
| `blrm_sim_one` | One cohort-by-cohort virtual trial | Inspect a single simulated path |
| `blrm_sim_n` | Repeated virtual trials | Generate trial records for operating characteristics |
| `blrm_sim_sum` | Long-format operating-characteristic summary | Summarize repeated simulations |
| `result_export` | Validated five-sheet XLSX export | Package analysis or simulation outputs |

## 5. Formal multi-chain interim analysis

Use `blrm_main_mc` for formal posterior inference.

```sas
%blrm_main_mc(
    datain=observed,
    dataout=blrm_posterior,
    ref_dose=12,
    mu1=-0.9,
    mu2=0,
    v1=4,
    v2=1,
    rho=0,
    alpha=0.05,
    nbi=2000,
    nmc=20000,
    nchain=4,
    simming=0,
    seed=9527,
    diagnostics=1,
    chain_diagout=blrm_chain,
    diagplots=0,
    ac_lags=%str(1 5 10 50),
    init_random=1,
    rhat_cutoff=1.01,
    postsumout=blrm_parameter_summary,
    diagout=blrm_rhat,
    debug=0
);
```

### Core parameters

| Parameter | Description |
|---|---|
| `datain` | Input DLT dataset |
| `dataout` | Pooled posterior draws from all chains |
| `ref_dose` | Positive reference dose |
| `mu1`, `mu2` | Prior means of `log_alpha` and `log_beta` |
| `v1`, `v2` | Prior variances |
| `rho` | Prior covariance |
| `alpha` | Tail probability for posterior intervals |
| `nbi` | Burn-in iterations per chain |
| `nmc` | Retained iterations per chain |
| `nchain` | Number of independent chains |
| `seed` | Seed for chain 1; subsequent chains use successive seeds |
| `diagnostics` | `1` saves within-chain diagnostic tables |
| `chain_diagout` | Prefix for pooled chain-level diagnostic datasets |
| `ac_lags` | Autocorrelation lags requested from `PROC MCMC` |
| `init_random` | `1` initializes chains independently from the prior |
| `rhat_cutoff` | Threshold used for the convergence flag |
| `postsumout` | Pooled posterior parameter summary |
| `diagout` | Classical Gelman-Rubin R-hat output |
| `debug` | `1` retains internal temporary datasets |

### Principal outputs

`dataout` contains the pooled posterior draws and identifies the originating
chain and within-chain iteration. `postsumout` contains pooled summaries for
`log_alpha` and `log_beta`. `diagout` contains classical Gelman-Rubin R-hat,
the selected cutoff, and a convergence flag.

When `diagnostics=1`, four pooled diagnostic datasets are also created from
the `chain_diagout` prefix:

- `<prefix>_ess`
- `<prefix>_mcse`
- `<prefix>_autocorr`
- `<prefix>_geweke`

Review convergence before using pooled draws for dose decisions. R-hat alone
is not sufficient: also inspect effective sample size, Monte Carlo standard
error, autocorrelation, Geweke results, and trace plots when needed.

## 6. Posterior toxicity and EWOC summaries

```sas
%blrm_stat(
    dose_list=%str(2,4,6,8,12),
    ref_dose=12,
    mcmcout=blrm_posterior,
    stat=%str(mean, std, q1, median, q3),
    dataout=posterior_toxicity,
    debug=0
);

%blrm_ewoc(
    dose_list=%str(2,4,6,8,12),
    ref_dose=12,
    mcmcout=blrm_posterior,
    dataout=ewoc_summary,
    ud=0.15,
    od=0.35,
    ewoc=0.25,
    debug=0
);
```

For each candidate dose, `blrm_ewoc` calculates:

- underdosing: `pi_d < ud`
- target toxicity: `ud <= pi_d < od`
- overdosing: `pi_d >= od`

A dose is EWOC-acceptable when its posterior overdosing probability is no
greater than `ewoc`. The dose with the largest target-toxicity probability is
a toxicity-based statistical recommendation, not an automatic RP2D. RP2D
selection requires cross-functional benefit-risk review.

## 7. Simulation workflow

Simulation requires two conceptually separate inputs:

1. an assumed true dose-toxicity scenario; and
2. the proposed operating design, including the prior, starting dose, cohort
   size, maximum sample size, EWOC threshold, and dose-skipping rule.

The truth dataset must contain `dose` and `true_tox`.

```sas
data truth;
    input dose true_tox;
    datalines;
1  0.10
5  0.15
10 0.24
12 0.32
15 0.63
;
run;
```

### One virtual trial

```sas
%blrm_sim_one(
    dose_list=%str(1,5,10,12,15),
    ref_dose=10,
    mu1=-0.9,
    mu2=0,
    v1=4,
    v2=1,
    rho=0,
    nbi=1000,
    nmc=10000,
    ud=0.16,
    od=0.33,
    ewoc=0.25,
    true_tox=truth,
    start_dose=10,
    max_sample=30,
    cohortn=3,
    skip=0,
    seed=9527,
    sim_batch=1,
    debug=0
);
```

Principal outputs are `DLT_RECORD`, `EWOC_RECORD`, and `SIM_RECORD`.

`skip=0` permits escalation by at most one dose level. `skip=1` permits the
implemented wider escalation step. The trial stops when the maximum sample
size is used or no candidate dose remains EWOC-acceptable.

### Repeated virtual trials

```sas
%blrm_sim_n(
    dose_list=%str(1,5,10,12,15),
    ref_dose=10,
    true_tox=truth,
    start_dose=10,
    max_sample=30,
    cohortn=3,
    skip=0,
    mu1=-0.9,
    mu2=0,
    v1=4,
    v2=1,
    rho=0,
    nbi=1000,
    nmc=10000,
    ud=0.16,
    od=0.33,
    ewoc=0.25,
    seed=9527,
    sim_time=1000,
    debug=0,
    time_out=time_summary,
    time_sum=1
);
```

Principal outputs are `DLT_ALL`, `EWOC_ALL`, `SIM_ALL`, and the requested
timing dataset. Simulation uses the single-chain engine within each virtual
trial; `blrm_main_mc` is reserved for formal inference rather than nested
inside repeated design simulation.

### Operating characteristics

```sas
%blrm_sim_sum(
    sim_data=sim_all,
    dlt_data=dlt_all,
    true_tox=truth,
    target_range_low=0.16,
    target_range_high=0.33,
    acceptable_low=0.05,
    out_summary=oc_summary,
    debug=0
);
```

The long-format output contains:

- `N_Simulations`
- `Prop_Correct_Dose`
- `Prop_Acceptable_Dose`
- `Prop_Over_Toxic_Dose`
- `Prop_No_Dose_Selected`
- `Avg_Sample_Used`
- `Avg_Patients_Over_Toxic`

All simulated trials remain in the denominator, including trials with no dose
selection. A trial without a DLT exposure record contributes zero patients to
the over-toxic exposure metric rather than being removed.

## 8. Excel export

```sas
%result_export(
    path=path/to/output,
    file_name=blrm_results,
    info=scenario_information,
    dlt_data=dlt_all,
    ewoc_data=ewoc_all,
    sim_data=sim_all,
    oc_data=oc_summary
);
```

The macro validates all five input datasets before creating the workbook and
writes `INFO`, `DLT`, `EWOC`, `SIM`, and `OC` sheets. An existing workbook with
the same name is replaced.

## 9. Independent R/JAGS comparison

`r/blrm_compare_rjags.R` implements the same likelihood, dose scaling, prior,
toxicity intervals, and EWOC rule directly in JAGS. The matched example uses
2 mg with 1/6 DLTs and 4 mg with 0/3 DLTs, then evaluates 2, 4, 6, 8, and 12 mg.

Install JAGS 4.x plus the R packages `rjags` and `coda`, then run from the
repository root:

```r
source("r/blrm_compare_rjags.R")
```

The script writes R/JAGS summaries, SAS-versus-R/JAGS differences, and MCMC
diagnostics to `output/`. Small numerical differences are expected because the
implementations use independent MCMC samplers.

## 10. Validation and operational use

The `tests/` directory contains eight SAS validation programs covering the
current macro interfaces, posterior draws and diagnostics, simulation record
aggregation, operating-characteristic denominators, and the five-sheet Excel
export.

Before study use:

- rerun the validation programs in the intended SAS environment;
- pre-specify the prior, toxicity intervals, EWOC threshold, dose-skipping
  convention, stopping rules, and simulation scenarios;
- review MCMC convergence before interpreting formal interim outputs; and
- perform independent programming and statistical validation required by the
  applicable study and organization.

This repository is provided for research and reproducibility. It is not a
validated production or regulatory system.

