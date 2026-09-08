# Data

The six CSV files reproduce the corresponding sheets in `extraction.xlsx`.

| File | One row represents |
|---|---|
| studies.csv | A study/population, including design and participant flow |
| reports.csv | A publication linked to a study, including secondary reports |
| conditions.csv | An intervention condition within a study/population |
| outcomes.csv | A reported or derived outcome summary or contrast |
| derivations.csv | A statistical derivation, its inputs, method, and assumptions |
| data_limitations.csv | A data limitation or source-verification finding and its analytical implication |

`study_id` links studies and reports; `population_id` distinguishes study strata; `condition_id` identifies an intervention condition; `contrast_id` groups related comparisons; `outcome_id` links outcome rows to derivations. These are relational keys, not participant identifiers. Derivation and limitation identifiers label scientific records.

In outcomes, `value` is the estimate described by `metric` and `result_form`. `dispersion_type` specifies what `dispersion_value` represents, including SD, SEM, SE, and paired-difference SD. `ci_low` and `ci_high` retain reported confidence limits. `p_value` is text so that inequalities such as `<0.001` remain explicit. `n` is the sample size associated with that result. `window_min` and `timepoint_min` are in minutes. `units_original`, `auc_definition`, `direction`, and `comparator` define the scale and comparison.

`data_status` distinguishes reported, derived, digitized, and author-supplied values. Provenance columns give the source location. Publication titles and DOIs are in `reports.csv`. Blank cells mean unavailable or not applicable, never zero. Unit uncertainties and graph digitization remain documented; source values are not silently standardized.

The [Enomoto source file](author_provided/Enomoto_2026/README.md) contains aggregate paired-test statistics used in four outcome rows.
