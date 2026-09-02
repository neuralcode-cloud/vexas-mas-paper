# Complement profiling

Produces Figure 2D and Figure S5.

**Code and data not yet added to this repository.**

## Expected layout

```
analysis/complement/run_complement.py   case vs controls
data/complement/complement.csv          long format, one row per subject x timepoint x analyte
```

`complement.csv` columns, enforced by `tests/test_data_contracts.py`: `subject`
(UPN*), `group`, `timepoint` (relative label, not a date), `analyte`, `value`,
`unit`, `scale` (`log10` or `linear`), `included`, `exclusion_reason`.

Estimators are importable from `analysis/cytokines_flow/singlecase.py`.
