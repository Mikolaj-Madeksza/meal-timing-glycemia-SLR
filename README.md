# Meal timing and glycemia: systematic review

Data, risk-of-bias assessments, and reproducible analyses for the systematic review of randomized meal-timing studies and acute glycemic responses.

- [extraction.xlsx](extraction.xlsx): study characteristics, publications, intervention conditions, outcomes, derivations, and data limitations.
- [risk_of_bias.xlsx](risk_of_bias.xlsx): detailed ROBUST-RCT assessments and a formula-linked summary.
- [meta_analysis](meta_analysis/README.md): analysis inputs, R code, tables, figures, and validation results.

The extraction contains 17 studies, 18 study/population records, 53 conditions, and 276 outcome records. Publication details and DOIs are in the Reports sheet. Source references identify the supporting table, figure, or passage. Study, population, condition, and outcome identifiers link related scientific records; they do not identify participants.

The only pooled comparison is Enomoto 2026 and Nakamura 2021: 180-minute CGM glucose iAUC after later versus earlier dinner. Other studies are reported individually. See the [analysis plan](meta_analysis/analysis_plan.md) for the estimand, methods, and limitations.

The Enomoto author-supplied file contains aggregate paired-test statistics only. It includes no individual participant data or participant identifiers.

Risk-of-bias counts and the highest core-item judgment per study are calculated from Assessments. The highest core-item judgment is a review-defined summary, not an overall judgment supplied by ROBUST-RCT. Optional item 7 is excluded from that calculation. Textual design considerations are explanatory annotations. See [risk-of-bias methods](risk_of_bias.md).

Software and documentation are provided under the existing [MIT license](LICENSE).

Eligibility correction (9 September 2026): Saad 2012 and Yadav 2023 were removed because only labelled-meal assessment order was randomized. Lu 2022 contributes the randomized no-preload lunch-delay comparison. Lu's outcome-window and capillary missing-data limitations are retained in the extraction and risk-of-bias evidence. The two-study delayed-dinner meta-analysis is unchanged.
