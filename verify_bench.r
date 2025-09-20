# ===== verify_bench_desktop.R =====

suppressPackageStartupMessages({
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Please install.packages('jsonlite')")
  if (!requireNamespace("coin", quietly = TRUE)) stop("Please install.packages('coin')")
})

# ----------------------------- JSONL loader -----------------------------

load_jsonl <- function(path) {
  # Open plain or gz automatically
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, warn = FALSE, encoding = "UTF-8")

  first_nonempty <- lines[trimws(lines) != ""][1]
  if (!is.na(first_nonempty) && startsWith(trimws(first_nonempty), "[")) {
    obj <- tryCatch(jsonlite::fromJSON(paste(lines, collapse = "\n"), simplifyVector = FALSE),
                    error = function(e) NULL)
    if (is.null(obj)) return(list())
    if (is.list(obj) && !is.null(obj$rows) && is.list(obj$rows)) return(obj$rows)
    if (is.list(obj) && !is.null(obj$records) && is.list(obj$records)) return(obj$records)
    if (is.list(obj)) return(obj)
    return(list())
  }

  # JSONL with possible log prefaces
  out <- list()
  for (s in lines) {
    s <- trimws(s)
    if (nchar(s) == 0) next
    if (substr(s, 1, 1) != "{") {
      i <- regexpr("\\{", s)
      if (i == -1) next
      s <- substring(s, i)
    }
    ok <- tryCatch({
      obj <- jsonlite::fromJSON(s, simplifyVector = FALSE)
      out[[length(out) + 1]] <- obj            
      TRUE
    }, error = function(e) FALSE)
    if (!ok) next
  }
  out
}

# ----------------------------- Stats helpers -----------------------------

fmt <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x) || is.nan(x) || is.infinite(x)) return("NA")
  formatC(as.numeric(x), digits = 3, format = "g")
}

paired_t <- function(y, x) {
  if (length(y) < 2) return(c(NA_real_, NA_real_))
  tt <- suppressWarnings(t.test(y, x, paired = TRUE))
  c(unname(tt$statistic), unname(tt$p.value))
}

# Wilcoxon signed-rank with Pratt zeros; W+ reported to mirror 'pratt' zero_method
wilcoxon_pratt <- function(y, x) {
  d <- y - x
  d <- d[!is.na(d)]
  if (length(d) == 0 || all(abs(d) < .Machine$double.eps)) return(list(W = 0, p = 1))
  absd <- abs(d)
  r_all <- rank(absd, ties.method = "average")
  Wplus <- sum(r_all[d > 0])

  nz <- d[d != 0]
  ties_present <- any(duplicated(abs(nz)))
  use_exact <- (!ties_present) && (length(nz) <= 25)  # SciPy-like 'auto' heuristic
  dist <- if (use_exact) "exact" else "asymptotic"

  wt <- tryCatch(coin::wilcoxsign_test(d ~ 1, zero.method = "Pratt", distribution = dist),
                 error = function(e) NULL)
  pval <- if (is.null(wt)) NA_real_ else as.numeric(coin::pvalue(wt))
  list(W = Wplus, p = pval)
}

cohen_dz <- function(d) {
  d <- d[!is.na(d)]
  if (length(d) < 2) return(NA_real_)
  sdv <- sd(d)
  if (!is.finite(sdv) || sdv == 0) return(NA_real_)
  mean(d) / sdv
}

rank_biserial_and_win <- function(d) {
  d <- d[!is.na(d)]
  if (length(d) == 0) return(c(NA_real_, NA_real_))
  n_pos <- sum(d > 0); n_neg <- sum(d < 0)
  denom <- n_pos + n_neg
  if (denom == 0) return(c(0, 0.5))
  r_rb <- (n_pos - n_neg) / denom
  win <- n_pos / denom
  c(r_rb, win)
}

bootstrap_ci_mean <- function(d, n_boot = 10000, seed = 0, alpha = 0.05) {
  d <- d[!is.na(d)]
  if (length(d) == 0) return(c(NA_real_, NA_real_))
  set.seed(seed)
  n <- length(d)
  boots <- replicate(n_boot, mean(d[sample.int(n, n, replace = TRUE)]))
  lo <- as.numeric(stats::quantile(boots, probs = alpha / 2, names = FALSE))
  hi <- as.numeric(stats::quantile(boots, probs = 1 - alpha / 2, names = FALSE))
  c(lo, hi)
}

# Holm step-down exactly like in the Python script:
holm_like_python <- function(keys, pvals) {
  valid_idx <- which(!is.na(pvals) & is.finite(pvals))
  m <- length(valid_idx)
  out <- rep(NA_real_, length(pvals))
  if (m == 0) { names(out) <- keys; return(out) }
  ord <- order(pvals[valid_idx], na.last = NA)
  ranks <- seq_along(ord)
  for (k in seq_along(ord)) {
    i <- valid_idx[ord[k]]
    padj <- min(1, (m - ranks[k] + 1) * pvals[i])
    out[i] <- padj
  }
  names(out) <- keys
  out
}

tier_from_bundle <- function(name) {
  s <- tolower(as.character(name))
  if (startsWith(s, "easy-")) return("easy")
  if (startsWith(s, "hard-")) return("hard")
  if (startsWith(s, "veryhard-")) return("veryhard")
  "unknown"
}

performance_profile <- function(tidy, metric_key, roles, step = 0.05) {
  thresholds <- round(seq(0, 1, by = step), 5)
  rows <- list()
  for (r in roles) {
    vals <- tidy[[metric_key]][tidy$role == r]
    vals <- vals[!is.na(vals)]
    n <- length(vals)
    for (tau in thresholds) {
      frac <- if (n > 0) mean(vals >= tau) else NA_real_
      rows[[length(rows)+1]] <- data.frame(role = r, threshold = tau, fraction = frac)
    }
  }
  do.call(rbind, rows)
}

# ----------------------------- Core analysis -----------------------------

analyze_once <- function(df, metric_key, alpha = 0.05, n_boot = 10000, seed = 0, systems = NULL) {
  roles <- sort(unique(df$role))
  if (!is.null(systems)) {
    want <- systems
    roles <- roles[roles %in% want]
    roles <- want[want %in% roles]  
  }
  if (length(roles) < 2) return("Not enough systems to compare.\n")

  wide <- reshape(df[, c("bundle", "role", metric_key)], timevar = "role", idvar = "bundle", direction = "wide")
  colnames(wide) <- sub(paste0("^", metric_key, "\\."), "", colnames(wide))
  rownames(wide) <- wide$bundle

  # Per-system summary
  lines <- c()
  lines <- c(lines, sprintf("# Bench Statistics (%s)", metric_key), "")
  lines <- c(lines, "## Per-system summary (over bundles)")
  lines <- c(lines, "| System | N | Mean | Median | Std | Min | Max |",
                    "|:--|--:|--:|--:|--:|--:|--:|")
  for (r in roles) {
    x <- as.numeric(wide[[r]])
    x <- x[!is.na(x)]
    if (length(x) == 0) {
      lines <- c(lines, sprintf("| %s | 0 |  |  |  |  |  |", r))
    } else {
      lines <- c(lines, sprintf("| %s | %d | %s | %s | %s | %s | %s |",
                                r, length(x), fmt(mean(x)), fmt(median(x)),
                                fmt(sd(x)), fmt(min(x)), fmt(max(x))))
    }
  }
  lines <- c(lines, "")

  # Latency summary if present
  if ("lat_total_ms" %in% names(df)) {
    lines <- c(lines, "## Per-system runtime (total ms)",
                      "| System | N | Mean | Median | Std | Min | Max |",
                      "|:--|--:|--:|--:|--:|--:|--:|")
    for (r in roles) {
      lat_r <- df[df$role == r & !is.na(df$lat_total_ms), c("bundle", "lat_total_ms")]
      merged <- merge(wide["bundle"], lat_r, by = "bundle", all.x = FALSE)
      x <- merged$lat_total_ms
      x <- x[!is.na(x)]
      if (length(x) == 0) {
        lines <- c(lines, sprintf("| %s | 0 |  |  |  |  |  |", r))
      } else {
        lines <- c(lines, sprintf("| %s | %d | %s | %s | %s | %s | %s |",
                                  r, length(x), fmt(mean(x)), fmt(median(x)),
                                  fmt(sd(x)), fmt(min(x)), fmt(max(x))))
      }
    }
    lines <- c(lines, "")
  }

  # Pairwise tests
  lines <- c(lines, "## Pairwise comparisons (paired, two-sided)")
  header <- "| A vs B | N | mean(B-A) | t | p_t | W | p_w | p_t_adj | p_w_adj | Cohen_dz | r_rb | win_rate(B>A) | 95% CI mean(B-A) |"
  lines <- c(lines, header, "|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|")

  pair_keys_t <- c(); pair_p_t <- c()
  pair_keys_w <- c(); pair_p_w <- c()
  pair_rows <- list()

  if (length(roles) >= 2) {
    for (i in seq_len(length(roles) - 1)) for (j in seq.int(i + 1, length(roles))) {
      a <- roles[i]; b <- roles[j]
      sub <- wide[!is.na(wide[[a]]) & !is.na(wide[[b]]), c("bundle", a, b)]
      n <- nrow(sub)
      key_t <- sprintf("%s__%s__t", a, b)
      key_w <- sprintf("%s__%s__w", a, b)
      if (n < 2) {
        pair_rows[[length(pair_rows)+1]] <- list(a=a,b=b,n=n,md=NA,t=NA,p_t=NA,W=NA,p_w=NA,
                                                 key_t=key_t,key_w=key_w,dz=NA,rrb=NA,win=NA,
                                                 lo=NA,hi=NA)
        next
      }
      x <- sub[[a]]; y <- sub[[b]]; d <- y - x
      md <- mean(d, na.rm = TRUE)
      pt <- paired_t(y, x)
      wres <- wilcoxon_pratt(y, x)
      dz <- cohen_dz(d)
      rb <- rank_biserial_and_win(d)
      ci <- bootstrap_ci_mean(d, n_boot = n_boot, seed = seed, alpha = alpha)

      pair_keys_t <- c(pair_keys_t, key_t); pair_p_t <- c(pair_p_t, pt[2])
      pair_keys_w <- c(pair_keys_w, key_w); pair_p_w <- c(pair_p_w, wres$p)

      pair_rows[[length(pair_rows)+1]] <- list(a=a,b=b,n=n,md=md,t=pt[1],p_t=pt[2],W=wres$W,p_w=wres$p,
                                               key_t=key_t,key_w=key_w,dz=dz,rrb=rb[1],win=rb[2],
                                               lo=ci[1],hi=ci[2])
    }
  }

  p_t_adj <- holm_like_python(pair_keys_t, pair_p_t)
  p_w_adj <- holm_like_python(pair_keys_w, pair_p_w)

  for (row in pair_rows) {
    lines <- c(lines, sprintf(
      "| %s vs %s | %d | %s | %s | %s | %s | %s | %s | %s | %s | %s | [%s,%s] |",
      row$a, row$b, row$n, fmt(row$md), fmt(row$t), fmt(row$p_t),
      fmt(row$W), fmt(row$p_w), fmt(p_t_adj[[row$key_t]]), fmt(p_w_adj[[row$key_w]]),
      fmt(row$dz), fmt(row$rrb), fmt(row$win), fmt(row$lo), fmt(row$hi)
    ))
  }

  lines <- c(lines, "",
             "Notes:",
             "- mean(B-A) > 0 means system B outperforms A on the chosen metric.",
             "- Cohen_dz is effect size for paired samples (mean diff / sd diff).",
             "- r_rb is the paired rank-biserial correlation; win_rate is fraction of bundles where B > A among non-ties.",
             "- p-values adjusted with Holm-Bonferroni per test family (t-tests, Wilcoxon) across all pairs.")
  paste(lines, collapse = "\n")
}

# ----------------------------- Public runner -----------------------------

verify_bench <- function(
  jsonl,
  metric = "f1_edges",
  systems = c("static","2R","4R"),
  alpha = 0.05,
  boot = 10000,
  seed = 0,
  by_tier = TRUE,
  out = NULL,                 
  csv = NULL,                 
  tidy_out = NULL,            
  delta_out = NULL,           
  pairs = c("static:2R","static:4R"),
  cdf_out = NULL,             
  lat_csv = NULL              
) {
  rows <- load_jsonl(jsonl)
  if (length(rows) == 0) stop("No JSON rows found in input.")

  # Flatten
  recs <- list()
  for (r in rows) {
    bundle <- r$bundle; role <- r$role; sc <- r$score
    if (is.null(bundle) || is.null(role)) next
    rec <- list(bundle = as.character(bundle), role = as.character(role))
    if (is.list(sc)) {
      for (k in names(sc)) {
        v <- sc[[k]]
        if (is.numeric(v) && !is.nan(v)) rec[[k]] <- as.numeric(v)
      }
    }
    lat <- r$latency
    if (is.list(lat)) {
      rec$lat_total_ms   <- if (!is.null(lat$total))   as.numeric(lat$total)   else NA_real_
      rec$lat_reader_ms  <- if (!is.null(lat$Reader))  as.numeric(lat$Reader)  else NA_real_
      rec$lat_mapper_ms  <- if (!is.null(lat$Mapper))  as.numeric(lat$Mapper)  else NA_real_
      rec$lat_writer_ms  <- if (!is.null(lat$Writer))  as.numeric(lat$Writer)  else NA_real_
      rec$lat_planner_ms <- if (!is.null(lat$Planner)) as.numeric(lat$Planner) else NA_real_
    }
    recs[[length(recs)+1]] <- rec
  }
  if (length(recs) == 0) stop("No usable records found (missing bundle/role/score).")
  df_all <- bind_records(recs)

  # Ensure types
  if (!"bundle" %in% names(df_all) || !"role" %in% names(df_all)) stop("Missing 'bundle' or 'role' columns.")
  df_all$bundle <- as.character(df_all$bundle)
  df_all$role   <- as.character(df_all$role)
  num_cols <- setdiff(names(df_all), c("bundle","role"))
  df_all[num_cols] <- lapply(df_all[num_cols], function(x) suppressWarnings(as.numeric(x)))

  metric_key <- metric
  if (!(metric_key %in% names(df_all))) {
    alt <- paste0("score.", metric)
    if (alt %in% names(df_all)) {
      metric_key <- alt
    } else {
      stop(sprintf("Metric '%s' not found. Available: %s",
                   metric, paste(sort(setdiff(names(df_all), c("bundle","role"))), collapse=", ")))
    }
  }

  if (!is.null(csv)) {
    wide <- reshape(df_all[, c("bundle", "role", metric_key)], timevar = "role", idvar = "bundle", direction = "wide")
    colnames(wide) <- sub(paste0("^", metric_key, "\\."), "", colnames(wide))
    wide <- wide[order(wide$bundle), ]
    dir.create(dirname(csv), showWarnings = FALSE, recursive = TRUE)
    write.csv(wide, csv, row.names = FALSE)
  }

  # Overall
  df_sub <- df_all[, c("bundle","role",metric_key, intersect(names(df_all), grep("^lat_", names(df_all), value = TRUE)))]
  report <- analyze_once(df_sub, metric_key, alpha, boot, seed, systems)

  # By tier
  if (isTRUE(by_tier)) {
    tiers <- c("easy","hard","veryhard")
    for (t in tiers) {
      mask <- startsWith(tolower(df_sub$bundle), paste0(t, "-"))
      if (any(mask)) {
        section <- analyze_once(df_sub[mask, ], metric_key, alpha, boot, seed, systems)
        report <- paste0(report, "\n\n---\n\n# ", tools::toTitleCase(t), " tier\n\n", section)
      }
    }
  }

  if (!is.null(out)) {
    dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
    writeLines(report, out, useBytes = TRUE)
    message(sprintf("Wrote %s", out))
  } else {
    cat(report, sep = "\n")
  }

  # Tidy CSV
  if (!is.null(tidy_out)) {
    tidy <- df_sub[, c("bundle","role",metric_key)]
    tidy$tier <- vapply(tidy$bundle, tier_from_bundle, character(1))
    colnames(tidy)[3] <- metric_key
    dir.create(dirname(tidy_out), showWarnings = FALSE, recursive = TRUE)
    write.csv(tidy, tidy_out, row.names = FALSE)
    message(sprintf("Wrote %s", tidy_out))
  }

  # Deltas
  if (!is.null(delta_out)) {
    pairs <- as.character(pairs)
    wide <- reshape(df_sub[, c("bundle","role",metric_key)], timevar = "role", idvar = "bundle", direction = "wide")
    colnames(wide) <- sub(paste0("^", metric_key, "\\."), "", colnames(wide))
    out_rows <- list()
    for (p in pairs) {
      if (!grepl(":", p, fixed = TRUE)) next
      sp <- strsplit(p, ":", fixed = TRUE)[[1]]
      a <- sp[1]; b <- sp[2]
      if (!(a %in% names(wide)) || !(b %in% names(wide))) next
      sub <- wide[!is.na(wide[[a]]) & !is.na(wide[[b]]), c("bundle", a, b)]
      if (nrow(sub) == 0) next
      tier <- vapply(sub$bundle, tier_from_bundle, character(1))
      diff <- sub[[b]] - sub[[a]]
      out_rows[[length(out_rows)+1]] <- data.frame(bundle = sub$bundle, pair = paste0(a, ":", b), diff = diff, tier = tier)
    }
    if (length(out_rows)) {
      deltas <- do.call(rbind, out_rows)
      dir.create(dirname(delta_out), showWarnings = FALSE, recursive = TRUE)
      write.csv(deltas, delta_out, row.names = FALSE)
      message(sprintf("Wrote %s", delta_out))
    }
  }

  # CDF profile
  if (!is.null(cdf_out)) {
    tidy <- df_sub[, c("bundle","role",metric_key)]
    roles <- unique(tidy$role)
    if (!is.null(systems)) {
      roles <- roles[roles %in% systems]
      roles <- systems[systems %in% roles]
    }
    cdf <- performance_profile(data.frame(tidy, check.names = FALSE), metric_key, roles)
    dir.create(dirname(cdf_out), showWarnings = FALSE, recursive = TRUE)
    write.csv(cdf, cdf_out, row.names = FALSE)
    message(sprintf("Wrote %s", cdf_out))
  }

  # Latency CSV
  if (!is.null(lat_csv)) {
    lat_cols <- grep("^lat_", names(df_all), value = TRUE)
    if (length(lat_cols)) {
      lat <- df_all[, c("bundle","role", lat_cols)]
      dir.create(dirname(lat_csv), showWarnings = FALSE, recursive = TRUE)
      write.csv(lat, lat_csv, row.names = FALSE)
      message(sprintf("Wrote %s", lat_csv))
    } else {
      message("No latency columns present in JSONL; skipped lat_csv.")
    }
  }

  invisible(report)
}

# ----------------------------- running the script -----------------------------
# jsonl_path <- file.choose()
# verify_bench(
#   jsonl = jsonl_path,
#   metric = "f1_edges",
#   systems = c("static","2R","4R"),
#   alpha = 0.05, boot = 10000, seed = 0, by_tier = TRUE,
#   out = "artifacts/bench_report_r.md",
#   csv = "artifacts/wide_f1_edges_r.csv",
#   tidy_out = "artifacts/tidy_f1_edges_r.csv",
#   delta_out = "artifacts/deltas_edges_r.csv", pairs = c("static:2R","static:4R"),
#   cdf_out = "artifacts/cdf_edges_r.csv",
#   lat_csv = "artifacts/latency_r.csv"
# )

# ===== end of file =====
