# BLRM in SAS for Single-Agent Phase I Dose Escalation

This repository provides a SAS implementation of the single-agent Bayesian logistic regression model (BLRM) for phase I dose-escalation studies.

## Repository Purpose

This repository is intended to provide:

- the SAS macro library for single-agent BLRM dose escalation
- a technical report describing the implementation
- a reusable starting point for posterior updating, EWOC-based decision support, and simulation

## Repository Structure

```text
blrm-sas-single-agent/
├─ README.md
├─ LICENSE
├─ CITATION.cff
├─ sas/
│  ├─ blrm_m2.sas
│  └─ blrm_export.sas
└─ report/
   └─ blrm_sas_technical_report.pdf
````

## Main Files

### `sas/blrm_m2.sas`

This is the main SAS macro library. It contains the core BLRM workflow, including functions for:

* posterior sampling with `PROC MCMC`
* posterior toxicity summarization
* EWOC-based dose recommendation
* single-trial execution
* repeated simulation
* simulation summarization

### `sas/blrm_export.sas`

This file contains export utilities used to write selected outputs to Excel-friendly formats.

### `report/blrm_sas_technical_report.pdf`

This report provides additional implementation details and usage examples beyond the repository overview.

## Software Requirements

* SAS 9.4 or later
* `PROC MCMC`

## Input Data Structure

The core BLRM workflow assumes a dose-level input dataset with the following variables:

* `dose`: numeric dose level
* `n`: number of treated subjects at that dose
* `dltn`: number of subjects with DLT at that dose

A simple example is:

```text
dose   n   dltn
1      3   0
5      6   1
```

## Minimal Usage Concept

A typical workflow is:

1. fit the BLRM model using the macro library in `sas/blrm_m2.sas`
2. summarize posterior toxicity probabilities for candidate doses
3. summarize underdosing, target toxicity, overdosing, and EWOC acceptability
4. if needed, run simulation to evaluate operating characteristics under prespecified truth scenarios

## Data Availability

This repository should not contain any restricted clinical trial data. Only toy examples or fully anonymized illustrative materials should be uploaded.

## Citation

If you use this repository, please cite the repository metadata provided in `CITATION.cff`.

## Contact

* Author: Zhanwen Peng
* Affiliation: AkesoBio
* Email: [316894075@qq.com](mailto:316894075@qq.com)

## License

This repository is released under the MIT License.

