# Proteomics — Olink Reveal

Produces Figure 2E.

## Inputs

In `../../data/proteomics/`:

| File | Contents |
|---|---|
| `olink_long.csv` | NPX table, 24 samples × 1033 assays, QC-passing rows only |
| `Reveal_Fixed_LOD.csv` | Olink fixed-LOD reference |
| `sample_metadata.csv` | De-identified sample sheet |

## Run

```bash
Rscript analysis/proteomics/run_pca.R
python -m pytest tests/test_proteomics.py -q
```

Requires R with `dplyr`, `tidyr`, `tibble`, `ggplot2`, `ggrepel`.

## Outputs

Written to `outputs/proteomics/`: `figure_2e_pca.png`, `pca_scores.csv`,
`pca_variance.csv`, `lod_per_assay.csv`, `assays_discarded.txt`,
`session_info.txt`.

## Notes

**LOD filtering.** The script prints the count under each of two thresholds:

| Rule | Discarded | Retained |
|---|---|---|
| below LOD in >95% of samples | 24 | 1009 |
| below LOD in 100% of samples | 15 | 1018 |

The manuscript reports 15 discarded and 1018 analysed. Both counts are asserted
in `tests/test_proteomics.py`.

**`original/PCA_Olink_data.R`** is the script as received from the study team,
unmodified apart from the removal of personal credit lines. `run_pca.R` differs
from it in five places:

1. Trailing comma removed from `return(list(..., discard_list = discardable, ))`.
2. `filter(Assay != unique(...))` replaced with `!Assay %in% ...`.
3. The placeholder input path replaced with `data/proteomics/olink_long.csv`.
4. `rownames()` assignment on a tibble dropped; labels come from
   `geom_text_repel`.
5. The unnamed row-index column in the CSV export is removed before use.

Verified against R 4.5.3.

**Sample metadata** columns: `sample_id`, `subject` (UPN*), `group`, `timepoint`
(T0/T1/T2), `age`, `sex`, `is_index_case`, `included`, `exclusion_reason`.
Derived from the study sample sheet using the timepoint logic of the original
script; all 24 sample identifiers match `olink_long.csv`. The source sheet
contains laboratory specimen accession numbers and is not in this repository.
