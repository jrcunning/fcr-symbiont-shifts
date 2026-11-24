# fcr-symbiont-shifts

Clean analysis of Florida coral symbiont composition shifts across sampling
eras, with sampling depth standardized before era contrasts.

## Project layout

- `code/analysis.Rmd`: main analysis notebook.
- `code/bleach2023_qpcr_exploration.Rmd`: exploratory analysis of the 2023 bleaching qPCR symbiont dataset.
- `R/symbiont_pipeline.R`: import, cleaning, depth-standardization, and test helpers.
- `metadata/`: sample metadata copied from the previous project.
- `data/raw/`: raw SymPortal output tables copied from the previous project.
- `data/processed/`: key reusable processed data products.
- `outputs/tables/`: analysis result tables.
- `outputs/figures/`: generated figures.

## Required files

Copy these files into the project before knitting `code/analysis.Rmd`:

- `metadata/decadalshifts_SST.xls`
- `data/raw/20240723T194127_Pawar_Run /post_med_seqs/500_20240730T092802_DBV_20240730T195050.seqs.absolute.abund_and_meta.txt`
- `data/raw/20250717T180136_conntr1/post_med_seqs/610_20250905T171203_DBV_20250910T163211.seqs.absolute.abund_and_meta.txt`
- `data/raw/experimental_shuffling_coef.csv`
- `data/raw/bleach2023SymQpcr.rds` for the separate qPCR exploration notebook.

The old Rmd also referenced type-profile tables, but the cleaned analysis does
not need them for the main clade-composition shift question.

Knitting the main analysis exports `data/processed/retained_sample_metadata.csv`
with one row for each sample retained after the subregion-by-species filter.
It includes all source metadata fields, derived host species codes and eras,
numeric latitude and longitude (decimal degrees), depth in feet, harmonized site
name, and post-MED read count. Collection dates use the source `YYYYMMDD`
format, and missing values are written as `NA`. Symbiont composition is not
included in this metadata export.

The retained sampling-location map also uses the R packages `sf`,
`rnaturalearth`, and `rnaturalearthhires` for the Florida outline. Knitting the
analysis writes a transparent, points-only overlay to
`outputs/figures/retained_sampling_locations_overlay.png`. Its extent is an
estimate of the reference image crop in Web Mercator; it may need adjustment to
align with a background whose exact extent or projection differs.
