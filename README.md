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
- `data/raw/20250702T202717_Conn_Run/post_med_seqs/606_20250710T083602_DBV_20250712T034645.seqs.absolute.abund_and_meta.txt`
- `data/raw/experimental_shuffling_coef.csv`
- `data/raw/bleach2023SymQpcr.rds` for the separate qPCR exploration notebook.

The old Rmd also referenced type-profile tables, but the cleaned analysis does
not need them for the main clade-composition shift question.
