# Analysis

The [analysis plan](analysis_plan.md) describes the primary comparison and sensitivity analyses.

Run these commands from the repository root using R (base packages only):

```sh
Rscript meta_analysis/assess_evidence.R
Rscript meta_analysis/summarize_studies.R
Rscript meta_analysis/run_meta_analysis.R
```

The scripts use the supplied CSV files and write:

- `outputs/evidence_assessment/`: effect availability, precision derivations, comparison compatibility, and uncertainty under assumed crossover correlations. These calculations do not pool studies or supply imputed variances to the primary model.
- `outputs/study_summary/`: individual study/population results with source units and provenance.
- `outputs/primary_meta_analysis/`: primary and sensitivity estimates, forest plots, and numerical REML validation.

The primary model uses Enomoto Day 2 and Nakamura Day 2. Enomoto Days 1, 3, and 4 and alternative Nakamura P-threshold assumptions are sensitivity analyses. Source values remain in their reported units; the pooled effects are converted to mmol/L × min. Positive values indicate greater glucose iAUC after later dinner.

To regenerate the CSV inputs from `extraction.xlsx` and check the workbooks and author-supplied statistics:

```sh
python3 -m pip install -r meta_analysis/requirements.txt
python3 meta_analysis/export_inputs.py
python3 meta_analysis/validate_data.py
```

The primary analysis requires R with Cairo graphics support and Times New Roman installed for the figures. Software versions from the latest run are recorded beside the outputs. See [data/README.md](data/README.md) for data structure and missing-value conventions.
