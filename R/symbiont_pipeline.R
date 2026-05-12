library(janitor)
library(readr)
library(readxl)
library(tidyverse)
library(vegan)

eras <- c("2001-03", "2019-21", "2024")
clades <- c("A", "B", "C", "D", "G")

required_files <- tribble(
  ~key, ~path, ~description,
  "metadata", "metadata/decadalshifts_SST.xls", "Sample metadata with sample_name, collection_date, collection_region, collection_subregion, collection_site, collection_depth_ft, host_genus, and host_species.",
  "pawar_variants", "data/raw/20240723T194127_Pawar_Run /post_med_seqs/500_20240730T092802_DBV_20240730T195050.seqs.absolute.abund_and_meta.txt", "SymPortal post-MED sequence absolute abundance table for the Pawar run.",
  "conn_variants", "data/raw/20250717T180136_conntr1/post_med_seqs/610_20250905T171203_DBV_20250910T163211.seqs.absolute.abund_and_meta.txt", "SymPortal post-MED sequence absolute abundance table for the Conn Dry Tortugas run.",
  "experimental_shuffling", "data/raw/experimental_shuffling_coef.csv", "Experimental shuffling coefficients and experimental clade-D summaries by coral species."
)

load_sample_metadata <- function(path) {
  readxl::read_xls(path) |>
    clean_names() |>
    mutate(
      collection_date = as.character(collection_date),
      year = readr::parse_integer(str_sub(collection_date, 1, 4)),
      era = case_when(
        year < 2005 ~ "2001-03",
        year >= 2023 ~ "2024",
        !is.na(year) ~ "2019-21",
        TRUE ~ NA_character_
      ),
      era = factor(era, levels = eras),
      depth_ft = parse_number(as.character(collection_depth_ft)),
      latitude = parse_number(as.character(collection_latitude)),
      longitude = parse_number(as.character(collection_longitude)),
      host_species_code = str_to_upper(str_c(str_sub(host_genus, 1, 1), str_sub(host_species, 1, 3)))
    )
}

read_symportal_variants <- function(path, run_id) {
  raw <- readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
  clean_names_vec <- janitor::make_clean_names(names(raw))
  sample_col <- match("sample_name", clean_names_vec)
  metadata_end_col <- match("collection_depth", clean_names_vec)

  if (is.na(sample_col) || is.na(metadata_end_col)) {
    stop("Could not find sample_name and collection_depth columns in: ", path)
  }

  metadata <- raw[, seq_len(metadata_end_col), drop = FALSE]
  names(metadata) <- clean_names_vec[seq_len(metadata_end_col)]
  metadata <- mutate(metadata, .row_id = row_number())

  abundance <- raw[, (metadata_end_col + 1):ncol(raw), drop = FALSE]
  names(abundance) <- names(raw)[(metadata_end_col + 1):ncol(raw)]

  metadata <- metadata |>
    filter(
      !is.na(sample_name),
      sample_name != "",
      !str_detect(str_to_lower(sample_name), "total")
    ) |>
    mutate(run_id = run_id)

  abundance <- abundance[metadata$.row_id, , drop = FALSE]
  abundance[] <- map(abundance, \(x) {
    if (is.character(x)) readr::parse_number(x) else as.numeric(x)
  })
  metadata <- select(metadata, -.row_id)

  list(
    samples = metadata,
    abundance = bind_cols(sample_name = metadata$sample_name, abundance) |>
      pivot_longer(-sample_name, names_to = "DIV", values_to = "count", values_drop_na = TRUE) |>
      mutate(run_id = run_id),
    taxonomy = tibble(
      DIV = names(abundance),
      clade = case_when(
        str_detect(DIV, "^[A-Z]") ~ str_extract(DIV, "^[A-Z]"),
        str_detect(DIV, "_[A-Z]$") ~ str_extract(DIV, "[A-Z]$"),
        TRUE ~ NA_character_
      )
    ) |>
      distinct()
  )
}

load_variant_runs <- function(run_files) {
  runs <- imap(run_files, ~ read_symportal_variants(.x, .y))

  list(
    samples = map_dfr(runs, "samples"),
    abundance = map_dfr(runs, "abundance"),
    taxonomy = map_dfr(runs, "taxonomy") |> distinct()
  )
}

load_experimental_shuffling <- function(path) {
  readr::read_csv(path, show_col_types = FALSE) |>
    clean_names() |>
    rename(
      species_code = species,
      experimental_n = n
    ) |>
    mutate(species_code = str_to_upper(species_code))
}

audit_sample_retention <- function(metadata, variants, min_reads = 100, clade_cols = clades) {
  sample_reads <- variants$abundance |>
    group_by(sample_name, run_id) |>
    summarize(post_med_reads = sum(count, na.rm = TRUE), .groups = "drop")

  clade_reads <- variants$abundance |>
    left_join(variants$taxonomy, by = "DIV") |>
    filter(clade %in% clade_cols) |>
    group_by(sample_name, run_id) |>
    summarize(target_clade_reads = sum(count, na.rm = TRUE), .groups = "drop")

  audited_samples <- variants$samples |>
    distinct(sample_name, run_id, .keep_all = TRUE) |>
    select(sample_name, run_id) |>
    left_join(sample_reads, by = c("sample_name", "run_id")) |>
    left_join(clade_reads, by = c("sample_name", "run_id")) |>
    left_join(
      metadata |>
        select(sample_name, collection_region, collection_subregion, era, host_species_code),
      by = "sample_name"
    ) |>
    mutate(
      post_med_reads = replace_na(post_med_reads, 0),
      target_clade_reads = replace_na(target_clade_reads, 0),
      retention_status = case_when(
        is.na(collection_region) ~ "dropped: metadata unmatched",
        post_med_reads < min_reads ~ "dropped: post-MED reads < threshold",
        collection_region != "FL" ~ "dropped: non-Florida sample",
        is.na(era) ~ "dropped: missing era",
        target_clade_reads <= 0 ~ "dropped: zero target clade reads",
        TRUE ~ "kept"
      )
    )

  list(
    stage_counts = audited_samples |>
      count(run_id, retention_status, collection_region, collection_subregion, era, name = "n") |>
      arrange(run_id, retention_status, collection_region, collection_subregion, era),
    dropped_samples = audited_samples |>
      filter(retention_status != "kept") |>
      arrange(run_id, retention_status, sample_name),
    duplicate_sample_names = audited_samples |>
      count(sample_name, name = "n_runs") |>
      filter(n_runs > 1) |>
      arrange(desc(n_runs), sample_name)
  )
}

prepare_sample_clade_composition <- function(metadata, variants, min_reads = 100) {
  sample_reads <- variants$abundance |>
    group_by(sample_name) |>
    summarize(post_med_reads = sum(count, na.rm = TRUE), .groups = "drop")

  variants$abundance |>
    left_join(variants$taxonomy, by = "DIV") |>
    filter(clade %in% clades) |>
    group_by(sample_name, clade) |>
    summarize(reads = sum(count, na.rm = TRUE), .groups = "drop") |>
    left_join(sample_reads, by = "sample_name") |>
    mutate(prop = reads / post_med_reads) |>
    select(sample_name, clade, prop) |>
    pivot_wider(names_from = clade, values_from = prop, values_fill = 0) |>
    left_join(sample_reads, by = "sample_name") |>
    left_join(
      metadata |>
        select(
          sample_name, host_genus, host_species, host_species_code,
          collection_region, collection_subregion, collection_site,
          collection_date, year, era, depth_ft, latitude, longitude
        ),
      by = "sample_name"
    ) |>
    filter(collection_region == "FL", post_med_reads >= min_reads, !is.na(era)) |>
    mutate(
      across(any_of(clades), ~ replace_na(.x, 0)),
      clade_total = rowSums(across(any_of(clades)))
    ) |>
    filter(clade_total > 0) |>
    mutate(across(any_of(clades), ~ .x / clade_total)) |>
    mutate(dominant_clade = names(pick(any_of(clades)))[max.col(as.matrix(pick(any_of(clades))), ties.method = "first")]) |>
    select(-clade_total)
}

make_depth_bins <- function(depth, n_bins = 6) {
  edges <- quantile(depth, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE)
  edges <- unique(as.numeric(edges))

  if (length(edges) < 3) {
    return(factor(rep(NA_character_, length(depth))))
  }

  cut(depth, breaks = edges, include.lowest = TRUE)
}

depth_weights <- function(group_df, target_bins) {
  target <- tibble(depth_bin = target_bins) |>
    filter(!is.na(depth_bin)) |>
    count(depth_bin, name = "n_target") |>
    mutate(p_target = n_target / sum(n_target))

  source <- group_df |>
    filter(!is.na(depth_bin)) |>
    count(depth_bin, name = "n_source") |>
    mutate(p_source = n_source / sum(n_source))

  weight_table <- target |>
    left_join(source, by = "depth_bin") |>
    mutate(
      p_source = replace_na(p_source, 0),
      weight = if_else(p_source > 0, p_target / p_source, 0)
    )

  group_df |>
    left_join(select(weight_table, depth_bin, weight), by = "depth_bin") |>
    mutate(weight = replace_na(weight, 0)) |>
    pull(weight)
}

weighted_composition <- function(group_df, weights = NULL, clade_cols = clades) {
  x <- as.matrix(select(group_df, all_of(clade_cols)))
  storage.mode(x) <- "double"

  if (is.null(weights)) {
    weights <- rep(1, nrow(group_df))
  }

  if (sum(weights, na.rm = TRUE) == 0) {
    weights <- rep(1, nrow(group_df))
  }

  weights <- weights / sum(weights, na.rm = TRUE)
  p <- colSums(x * weights, na.rm = TRUE)

  if (sum(p) == 0) {
    set_names(rep(NA_real_, length(clade_cols)), clade_cols)
  } else {
    set_names(p / sum(p), clade_cols)
  }
}

hellinger_distance <- function(p, q) {
  sqrt(0.5 * sum((sqrt(p) - sqrt(q))^2))
}

total_variation_distance <- function(p, q) {
  0.5 * sum(abs(p - q))
}

jensen_shannon_distance <- function(p, q) {
  midpoint <- 0.5 * (p + q)
  kl <- function(a, b) {
    keep <- a > 0
    sum(a[keep] * (log2(a[keep]) - log2(b[keep])))
  }
  sqrt(0.5 * kl(p, midpoint) + 0.5 * kl(q, midpoint))
}

correlation_summary <- function(x, y, cor_name, p_name) {
  complete <- complete.cases(x, y)
  x <- x[complete]
  y <- y[complete]

  if (length(x) < 3 || sd(x) == 0 || sd(y) == 0) {
    return(tibble(n_species = length(x), !!cor_name := NA_real_, !!p_name := NA_real_))
  }

  test <- cor.test(x, y)
  tibble(n_species = length(x), !!cor_name := unname(test$estimate), !!p_name := test$p.value)
}

distance_correlation_summary <- function(df) {
  hellinger <- correlation_summary(
    df$Hellinger,
    df$mean_shuff_coef,
    "cor_hellinger_shuff_coef",
    "p_hellinger_shuff_coef"
  )
  jsd <- correlation_summary(
    df$JSD,
    df$mean_shuff_coef,
    "cor_jsd_shuff_coef",
    "p_jsd_shuff_coef"
  ) |>
    select(-n_species)

  bind_cols(hellinger, jsd)
}

standardized_era_composition <- function(data, group_cols = c("collection_subregion", "host_species_code"),
                                         min_n = 4, n_bins = 6, clade_cols = clades,
                                         no_depth_standardization_subregions = "Dry Tortugas") {
  data |>
    group_by(across(all_of(group_cols))) |>
    group_modify(\(.x, .y) {
      df <- .x
      subregion <- if ("collection_subregion" %in% names(.y)) .y$collection_subregion else NA_character_
      skip_depth_standardization <- !is.na(subregion) && subregion %in% no_depth_standardization_subregions
      has_complete_depth <- all(!is.na(df$depth_ft))
      depth_bins <- if (has_complete_depth) make_depth_bins(df$depth_ft, n_bins = n_bins) else factor(rep(NA_character_, nrow(df)))
      use_depth_standardization <- !skip_depth_standardization && has_complete_depth && !all(is.na(depth_bins))

      if (use_depth_standardization) {
        df <- mutate(df, depth_bin = depth_bins)
        target_bins <- df$depth_bin
        standardization_method <- "Depth-standardized"
      } else if (skip_depth_standardization) {
        standardization_method <- "Not depth-standardized: Dry Tortugas"
      } else {
        standardization_method <- "Not depth-standardized"
      }

      map_dfr(eras, \(era_i) {
        era_df <- filter(df, era == era_i)
        n_samples <- nrow(era_df)

        if (n_samples < min_n) {
          comp <- set_names(rep(NA_real_, length(clade_cols)), clade_cols)
        } else if (use_depth_standardization) {
          comp <- weighted_composition(era_df, depth_weights(era_df, target_bins), clade_cols)
        } else {
          comp <- weighted_composition(era_df, clade_cols = clade_cols)
        }

        tibble(
          era = era_i,
          n = n_samples,
          depth_standardized = use_depth_standardization,
          standardization_method = standardization_method,
          !!!as.list(comp)
        )
      })
    }) |>
    ungroup() |>
    mutate(era = factor(era, levels = eras))
}

compare_eras <- function(composition, clade_cols = clades) {
  contrasts <- tribble(
    ~era1, ~era2,
    "2001-03", "2019-21",
    "2019-21", "2024",
    "2001-03", "2024"
  )

  composition |>
    group_by(collection_subregion, host_species_code) |>
    group_modify(\(.x, .y) {
      map_dfr(seq_len(nrow(contrasts)), \(i) {
        e1 <- contrasts$era1[i]
        e2 <- contrasts$era2[i]
        p1 <- filter(.x, era == e1)
        p2 <- filter(.x, era == e2)
        n1 <- if (nrow(p1) == 1) p1$n else NA_integer_
        n2 <- if (nrow(p2) == 1) p2$n else NA_integer_
        depth_standardized <- if (nrow(p1) == 1) p1$depth_standardized else NA
        standardization_method <- if (nrow(p1) == 1) p1$standardization_method else NA_character_

        if (nrow(p1) != 1 || nrow(p2) != 1 || anyNA(select(p1, all_of(clade_cols))) || anyNA(select(p2, all_of(clade_cols)))) {
          return(tibble(
            era1 = e1,
            era2 = e2,
            n1 = n1,
            n2 = n2,
            depth_standardized = depth_standardized,
            standardization_method = standardization_method,
            TV = NA_real_,
            Hellinger = NA_real_,
            JSD = NA_real_
          ))
        }

        v1 <- as.numeric(select(p1, all_of(clade_cols)))
        v2 <- as.numeric(select(p2, all_of(clade_cols)))

        tibble(
          era1 = e1,
          era2 = e2,
          n1 = n1,
          n2 = n2,
          depth_standardized = depth_standardized,
          standardization_method = standardization_method,
          TV = total_variation_distance(v1, v2),
          Hellinger = hellinger_distance(v1, v2),
          JSD = jensen_shannon_distance(v1, v2)
        )
      })
    }) |>
    ungroup()
}

compare_its2_to_experimental <- function(composition, distances, experimental) {
  composition_2024 <- composition |>
    filter(era == "2024") |>
    select(
      collection_subregion,
      species_code = host_species_code,
      n_its2_2024 = n,
      its2_prop_d_2024 = D
    )

  distance_joined <- distances |>
    mutate(
      species_code = host_species_code,
      contrast = paste(era1, "to", era2)
    ) |>
    left_join(experimental, by = "species_code")

  composition_joined <- composition_2024 |>
    left_join(experimental, by = "species_code")

  distance_correlations <- distance_joined |>
    filter(!is.na(mean_shuff_coef), !is.na(Hellinger)) |>
    group_by(collection_subregion, contrast) |>
    group_modify(\(.x, .y) distance_correlation_summary(.x)) |>
    ungroup()

  composition_correlations <- composition_joined |>
    filter(!is.na(mean_prop_d), !is.na(its2_prop_d_2024)) |>
    group_by(collection_subregion) |>
    group_modify(\(.x, .y) {
      correlation_summary(
        .x$its2_prop_d_2024,
        .x$mean_prop_d,
        "cor_its2_d_experimental_d",
        "p_its2_d_experimental_d"
      )
    }) |>
    ungroup()

  list(
    distance = distance_joined,
    composition = composition_joined,
    distance_correlations = distance_correlations,
    composition_correlations = composition_correlations
  )
}

combine_regions_by_species <- function(composition, regions = c("Upper Keys", "Dry Tortugas"),
                                       clade_cols = clades, region_weighting = c("equal", "sample_n")) {
  region_weighting <- match.arg(region_weighting)

  composition |>
    filter(collection_subregion %in% regions) |>
    filter(if_any(all_of(clade_cols), ~ !is.na(.x))) |>
    group_by(host_species_code, era) |>
    group_modify(\(.x, .y) {
      if (region_weighting == "equal") {
        weights <- rep(1, nrow(.x))
      } else {
        weights <- .x$n
      }

      weights <- weights / sum(weights, na.rm = TRUE)
      x <- as.matrix(select(.x, all_of(clade_cols)))
      combined <- colSums(x * weights, na.rm = TRUE)
      combined <- combined / sum(combined)

      tibble(
        n_total = sum(.x$n, na.rm = TRUE),
        n_regions = n_distinct(.x$collection_subregion),
        regions_included = str_c(sort(unique(.x$collection_subregion)), collapse = "; "),
        region_weighting = region_weighting,
        standardization_context = str_c(sort(unique(.x$standardization_method)), collapse = "; "),
        !!!as.list(set_names(combined, clade_cols))
      )
    }) |>
    ungroup() |>
    mutate(era = factor(era, levels = eras))
}

compare_eras_by_species <- function(composition, clade_cols = clades) {
  contrasts <- tribble(
    ~era1, ~era2,
    "2001-03", "2019-21",
    "2019-21", "2024",
    "2001-03", "2024"
  )

  composition |>
    group_by(host_species_code) |>
    group_modify(\(.x, .y) {
      map_dfr(seq_len(nrow(contrasts)), \(i) {
        e1 <- contrasts$era1[i]
        e2 <- contrasts$era2[i]
        p1 <- filter(.x, era == e1)
        p2 <- filter(.x, era == e2)
        n1 <- if (nrow(p1) == 1) p1$n_total else NA_integer_
        n2 <- if (nrow(p2) == 1) p2$n_total else NA_integer_
        regions_era1 <- if (nrow(p1) == 1) p1$regions_included else NA_character_
        regions_era2 <- if (nrow(p2) == 1) p2$regions_included else NA_character_

        if (nrow(p1) != 1 || nrow(p2) != 1 || anyNA(select(p1, all_of(clade_cols))) || anyNA(select(p2, all_of(clade_cols)))) {
          return(tibble(
            era1 = e1,
            era2 = e2,
            n1 = n1,
            n2 = n2,
            regions_era1 = regions_era1,
            regions_era2 = regions_era2,
            TV = NA_real_,
            Hellinger = NA_real_,
            JSD = NA_real_
          ))
        }

        v1 <- as.numeric(select(p1, all_of(clade_cols)))
        v2 <- as.numeric(select(p2, all_of(clade_cols)))

        tibble(
          era1 = e1,
          era2 = e2,
          n1 = n1,
          n2 = n2,
          regions_era1 = regions_era1,
          regions_era2 = regions_era2,
          TV = total_variation_distance(v1, v2),
          Hellinger = hellinger_distance(v1, v2),
          JSD = jensen_shannon_distance(v1, v2)
        )
      })
    }) |>
    ungroup()
}

compare_across_region_to_experimental <- function(composition, distances, experimental) {
  composition_2024 <- composition |>
    filter(era == "2024") |>
    select(
      species_code = host_species_code,
      n_its2_2024 = n_total,
      n_regions_2024 = n_regions,
      regions_2024 = regions_included,
      its2_prop_d_2024 = D
    )

  distance_joined <- distances |>
    mutate(
      species_code = host_species_code,
      contrast = paste(era1, "to", era2)
    ) |>
    left_join(experimental, by = "species_code")

  composition_joined <- composition_2024 |>
    left_join(experimental, by = "species_code")

  distance_correlations <- distance_joined |>
    filter(!is.na(mean_shuff_coef), !is.na(Hellinger)) |>
    group_by(contrast) |>
    group_modify(\(.x, .y) distance_correlation_summary(.x)) |>
    ungroup()

  composition_correlations <- composition_joined |>
    filter(!is.na(mean_prop_d), !is.na(its2_prop_d_2024))

  composition_correlations <- correlation_summary(
    composition_correlations$its2_prop_d_2024,
    composition_correlations$mean_prop_d,
    "cor_its2_d_experimental_d",
    "p_its2_d_experimental_d"
  )

  list(
    distance = distance_joined,
    composition = composition_joined,
    distance_correlations = distance_correlations,
    composition_correlations = composition_correlations
  )
}

partial_pool_all_region_composition <- function(sample_data,
                                                regions = c("Upper Keys", "Dry Tortugas", "Marquesas"),
                                                clade_cols = clades) {
  model_data <- sample_data |>
    filter(collection_subregion %in% regions, era %in% eras) |>
    mutate(
      host_species_code = factor(host_species_code),
      collection_subregion = factor(collection_subregion, levels = regions),
      era = factor(era, levels = eras),
      era_num = as.numeric(era) - 1
    )

  observed_species_regions <- model_data |>
    distinct(host_species_code, collection_subregion)

  prediction_grid <- observed_species_regions |>
    tidyr::expand_grid(era = factor(eras, levels = eras)) |>
    mutate(
      host_species_code = factor(host_species_code, levels = levels(model_data$host_species_code)),
      collection_subregion = factor(collection_subregion, levels = levels(model_data$collection_subregion)),
      era_num = as.numeric(era) - 1
    )

  predict_one_clade <- function(clade) {
    dat <- model_data |>
      mutate(
        clade_reads = pmin(post_med_reads, pmax(0, round(.data[[clade]] * post_med_reads))),
        other_reads = post_med_reads - clade_reads
      )

    fit <- mgcv::gam(
      cbind(clade_reads, other_reads) ~
        era_num +
        collection_subregion +
        s(host_species_code, bs = "re") +
        s(host_species_code, by = era_num, bs = "re"),
      data = dat,
      family = quasibinomial(),
      method = "REML"
    )

    prediction_grid |>
      mutate(!!clade := as.numeric(predict(fit, newdata = prediction_grid, type = "response"))) |>
      select(host_species_code, collection_subregion, era, all_of(clade))
  }

  predicted_region_composition <- map(clade_cols, predict_one_clade) |>
    reduce(left_join, by = c("host_species_code", "collection_subregion", "era")) |>
    mutate(
      across(all_of(clade_cols), ~ replace_na(.x, 0)),
      pred_sum = rowSums(across(all_of(clade_cols)))
    ) |>
    filter(pred_sum > 0) |>
    mutate(across(all_of(clade_cols), ~ .x / pred_sum)) |>
    select(-pred_sum)

  observed_counts <- model_data |>
    count(host_species_code, era, name = "n_total")

  predicted_region_composition |>
    group_by(host_species_code, era) |>
    summarize(
      n_regions = n_distinct(collection_subregion),
      regions_included = str_c(sort(unique(as.character(collection_subregion))), collapse = "; "),
      region_weighting = "partial pooling; equal observed-region prediction weights",
      standardization_context = "GAM one-vs-rest clade models with species random intercepts/slopes",
      across(all_of(clade_cols), mean),
      .groups = "drop"
    ) |>
    left_join(observed_counts, by = c("host_species_code", "era")) |>
    mutate(
      n_total = replace_na(n_total, 0L),
      era = factor(era, levels = eras)
    ) |>
    relocate(n_total, .after = era)
}

permanova_by_group <- function(data, group_cols = c("collection_subregion", "host_species_code"),
                               min_n = 4, clade_cols = clades, permutations = 999,
                               no_depth_adjustment_subregions = "Dry Tortugas") {
  data |>
    group_by(across(all_of(group_cols))) |>
    group_modify(\(.x, .y) {
      df <- .x
      subregion <- if ("collection_subregion" %in% names(.y)) .y$collection_subregion else NA_character_
      skip_depth_adjustment <- !is.na(subregion) && subregion %in% no_depth_adjustment_subregions
      has_complete_depth <- all(!is.na(df$depth_ft))
      depth_bins <- if (has_complete_depth) make_depth_bins(df$depth_ft) else factor(rep(NA_character_, nrow(df)))
      use_depth_adjustment <- !skip_depth_adjustment && has_complete_depth && !all(is.na(depth_bins)) && n_distinct(depth_bins) > 1

      if (use_depth_adjustment) {
        df <- mutate(df, depth_bin = depth_bins)
      }

      era_counts <- df |>
        count(era) |>
        filter(n >= min_n)

      if (nrow(era_counts) < 2) {
        return(tibble(
          n = nrow(df),
          eras_tested = nrow(era_counts),
          depth_adjusted = use_depth_adjustment,
          test_formula = if_else(use_depth_adjustment, "composition ~ era + depth_bin", "composition ~ era"),
          r2_era = NA_real_,
          p_era = NA_real_,
          r2_depth_bin = NA_real_,
          p_depth_bin = NA_real_
        ))
      }

      df <- filter(df, era %in% era_counts$era)
      if (use_depth_adjustment) {
        fit <- vegan::adonis2(
          as.matrix(select(df, all_of(clade_cols))) ~ era + depth_bin,
          data = df,
          method = "bray",
          permutations = permutations
        )

        r2_depth_bin <- fit$R2[2]
        p_depth_bin <- fit$`Pr(>F)`[2]
      } else {
        fit <- vegan::adonis2(
          as.matrix(select(df, all_of(clade_cols))) ~ era,
          data = df,
          method = "bray",
          permutations = permutations
        )

        r2_depth_bin <- NA_real_
        p_depth_bin <- NA_real_
      }

      tibble(
        n = nrow(df),
        eras_tested = n_distinct(df$era),
        depth_adjusted = use_depth_adjustment,
        test_formula = if_else(use_depth_adjustment, "composition ~ era + depth_bin", "composition ~ era"),
        r2_era = fit$R2[1],
        p_era = fit$`Pr(>F)`[1],
        r2_depth_bin = r2_depth_bin,
        p_depth_bin = p_depth_bin
      )
    }) |>
    ungroup()
}
