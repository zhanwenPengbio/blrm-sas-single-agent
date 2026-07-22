# BLRM in SAS for Single-Agent Phase I Dose Escalation

This repository provides a SAS implementation of the single-agent Bayesian logistic regression model (BLRM) for phase I oncology dose-escalation studies. It supports posterior updating, EWOC-based dose decision support, independent multi-chain analysis with MCMC diagnostics, and simulation-based evaluation of operating characteristics.

## Scope

The code is intended to support statistical analysis and simulation. A dose with the highest posterior probability of lying in the target toxicity interval is not automatically the recommended phase II dose (RP2D). RP2D selection should integrate safety, efficacy, pharmacokinetics, pharmacodynamics, dose intensity, and clinical judgment through cross-functional review.

## Model

For dose \(d\) and reference dose \(d^*\), the model is

```text
logit(pi_d) = log_alpha + exp(log_beta) * log(d / d*)
```

with the prior

```text
(log_alpha, log_beta)' ~ MVN(mu, Sigma)

Sigma = [ v1   rho ]
        [ rho  v2  ]
```

In this implementation, `v1` and `v2` are prior variances and `rho` is the prior covariance, not a correlation coefficient. The supplied covariance matrix must be positive definite.

## Repository Structure

```text
blrm-sas-single-agent/
├── README.md
├── LICENSE
├── CITATION.cff
├── docs/
│   └── user-guide.md
├── r/
│   └── blrm_compare_rjags.R
├── sas/
│   ├── blrm_m2.sas
│   └── blrm_export.sas
├── tests/
│   └── test_01_... through test_08_...
└── report/
    └── blrm_sas_technical_report.pdf
```

See the [BLRM SAS Macro User Guide](docs/user-guide.md) for the complete
interim-analysis, diagnostic, simulation, and export workflow.

## Main SAS Macros

| Macro | Purpose |
|---|---|
| `blrm_main` | Runs a single `PROC MCMC` chain and optionally saves ESS, MCSE, autocorrelation, and Geweke diagnostics. |
| `blrm_main_mc` | Runs independent chains, pools posterior draws, saves chain-level diagnostics, and calculates classical Gelman-Rubin R-hat for `log_alpha` and `log_beta`. |
| `blrm_stat` | Summarizes posterior toxicity probabilities at specified dose levels. |
| `blrm_ewoc` | Calculates posterior underdosing, target-toxicity, and overdosing probabilities and applies the EWOC rule. |
| `blrm_sim_one` | Simulates one cohort-by-cohort BLRM trial. |
| `blrm_sim_n` | Repeats `blrm_sim_one` and combines trial-level records. |
| `blrm_sim_sum` | Summarizes dose-selection and patient-exposure operating characteristics. |
| `result_export` | Validates and exports the current `INFO`, `DLT`, `EWOC`, `SIM`, and long-format `OC` data sets to one XLSX workbook. |

The internal helper macros `make_dose_map` and `make_dose_truth_map` preserve the order of the supplied dose list and align simulated toxicity probabilities by numeric dose value.

## Software Requirements

- SAS 9.4 or later
- `PROC MCMC`
- ODS diagnostic tables for ESS, MCSE, autocorrelation, and Geweke output when these diagnostics are requested

## Input Data

The fitting macros expect one record per evaluated dose level:

| Variable | Description |
|---|---|
| `dose` | Numeric dose level |
| `n` | Number of treated participants |
| `dltn` | Number of participants with a dose-limiting toxicity |

Example:

```sas
data observed;
    input dose n dltn;
    datalines;
1 3 0
5 6 1
;
run;
```

## Formal Multi-Chain Analysis

Use `blrm_main_mc` for formal posterior inference and convergence assessment.

```sas
%include "path/to/sas/blrm_m2.sas";

%blrm_main_mc(
    datain=observed,
    dataout=blrm_posterior,
    ref_dose=5,
    mu1=-1.3863,
    mu2=0,
    v1=4,
    v2=1,
    rho=0,
    nbi=2000,
    nmc=20000,
    nchain=4,
    diagnostics=1,
    seed=9527,
    rhat_cutoff=1.01,
    postsumout=blrm_parameter_summary,
    diagout=blrm_rhat
);
```

The numerical prior values above are illustrative only and are not design recommendations.

Principal outputs include:

- `blrm_posterior`: pooled posterior draws with chain identifiers
- `blrm_parameter_summary`: pooled summaries for `log_alpha` and `log_beta`
- `blrm_rhat`: classical Gelman-Rubin R-hat and convergence flags
- `blrm_posterior_chain_diag_ess`
- `blrm_posterior_chain_diag_mcse`
- `blrm_posterior_chain_diag_autocorr`
- `blrm_posterior_chain_diag_geweke`

When `init_random=1`, each chain is initialized independently from the prior. Chain seeds are generated from the supplied base seed.

## Posterior Toxicity and EWOC Summaries

```sas
%blrm_stat(
    dose_list=%str(1,2,5,10),
    ref_dose=5,
    mcmcout=blrm_posterior,
    dataout=posterior_toxicity
);

%blrm_ewoc(
    dose_list=%str(1,2,5,10),
    ref_dose=5,
    mcmcout=blrm_posterior,
    dataout=ewoc_summary,
    ud=0.16,
    od=0.33,
    ewoc=0.25
);
```

For each dose, `blrm_ewoc` reports posterior probabilities for:

- underdosing: `pi_d < ud`
- target toxicity: `ud <= pi_d < od`
- overdosing: `pi_d >= od`

A dose satisfies the implemented EWOC rule when its posterior overdosing probability is no greater than `ewoc`.

## Simulation Workflow

Simulation can be used to:

1. predict trial behavior under toxicity scenarios supplied by the study team; and
2. compare the operating characteristics of BLRM with alternative dose-escalation designs.

The expected workflow is:

```text
true toxicity scenario
        |
        v
blrm_sim_one -> cohort data and EWOC decisions for one trial
        |
        v
blrm_sim_n   -> repeated trial records
        |
        v
blrm_sim_sum -> operating-characteristic summary
```

`blrm_sim_sum` returns the following metrics in long format:

- `N_Simulations`
- `Prop_Correct_Dose`
- `Prop_Acceptable_Dose`
- `Prop_Over_Toxic_Dose`
- `Prop_No_Dose_Selected`
- `Avg_Sample_Used`
- `Avg_Patients_Over_Toxic`

The summary keeps all simulated trials in the denominator, including trials for which no dose is selected. It also validates the toxicity-scenario mapping before calculating dose-selection and exposure metrics.

## Independent R/JAGS Check

The script [`r/blrm_compare_rjags.R`](r/blrm_compare_rjags.R) implements the
same BLRM directly in JAGS for an external implementation check. It reproduces
the matched comparison scenario used in the project paper and presentation,
calculates posterior underdosing, target-toxicity, and overdosing probabilities,
applies the static EWOC rule, and reports R-hat and effective sample size.

Requirements are JAGS 4.x and the R packages `rjags` and `coda`. The script
runs eight independent chains and writes pooled region probabilities, chain-level
stability summaries, MCMC diagnostics, and chain-specific results as four CSV files
in the current R working directory.

## Excel Export

Run `sas/blrm_export.sas` after producing the analysis and simulation outputs.

```sas
%include "path/to/sas/blrm_export.sas";

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

The macro validates all five input data sets before writing the workbook. It
creates `INFO`, `DLT`, `EWOC`, `SIM`, and `OC` sheets. If a workbook with the
same name already exists, it is replaced.

## Recommended Use

- Use `blrm_main` for simple runs, simulation internals, and detailed inspection of a single chain.
- Use `blrm_main_mc` for formal analysis requiring multiple independent chains and convergence diagnostics.
- Set `debug=1` only when intermediate data sets are needed for development or validation.
- Pre-specify the prior, toxicity intervals, EWOC threshold, dose-skipping rule, and simulation scenarios before interpreting results.

## Data Availability

This repository should contain no restricted clinical trial data. Use toy examples, simulated data, or fully anonymized illustrative inputs only.

## Disclaimer

This implementation is provided for research and reproducibility. It should be independently validated for the intended computing environment and trial protocol before operational use.
