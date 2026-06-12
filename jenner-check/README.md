# Jenner compatibility bundles

This directory was added by a pull request from the
[Jenner](https://jenneranalytics.com) project. Each `tNNN_*` subdirectory is a
small, self-contained SAS bundle built from code in this repository, so you can
see your own SAS running on Jenner.

## What's in here

```
jenner-check/
├── README.md            # this file
├── run_jenner.sh        # mac / linux runner (curl)
├── run_jenner.bat       # windows runner
├── run_jenner.sas       # run from base SAS (PROC HTTP, SAS 9.4 M5+)
└── tNNN_<name>/
    ├── script.sas       # the SAS under test, adapted from this repo
    ├── autoexec.sas     # options + the sample data the script reads
    ├── expected.json    # the fields pinned from a passing run
    └── expected/        # human-readable snapshot of that run
        ├── log.txt
        ├── output.txt
        └── files.md     # links to any generated files / data sets
```

## How to run it

From inside `jenner-check/`, against the hosted API:

```bash
./run_jenner.sh --all          # run every bundle
./run_jenner.sh t001_result_export   # run just one
```

On Windows use `run_jenner.bat`; from base SAS (9.4 M5+):

```sas
%include 'run_jenner.sas';
%jenner_check_all();
```

Each bundle posts `autoexec.sas` + `script.sas` to
`https://api.jenneranalytics.com/v1/run` and prints the status, log and
listing. You can also paste any `script.sas` into the collaborative workspace
at [jenneranalytics.com](https://jenneranalytics.com) and run it there.

## The bundles

| bundle | what it runs |
|--------|--------------|
| `t001_result_export` | the `%result_export` utility from `sas/blrm_export.sas`, writing the simulation outputs to an Excel workbook |
| `t002_posterior_toxicity` | the BLRM posterior dose-toxicity transform from `sas/blrm_m2.sas`, summarized per dose with PROC MEANS |

## Don't want future PRs from us?

Reply with `no-more-prs` anywhere in a comment, or open an issue titled
`jenner-check: opt out`, and we'll stop.
