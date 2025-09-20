# Bench Statistics (f1_edges)

## Per-system summary (over bundles)

| System |   N |  Mean | Median |   Std |   Min | Max |
| :----- | --: | ----: | -----: | ----: | ----: | --: |
| static |  90 | 0.684 |    0.5 | 0.246 | 0.286 |   1 |
| 2R     |  90 | 0.936 |      1 | 0.131 |   0.5 |   1 |
| 4R     |  90 | 0.953 |      1 | 0.111 | 0.667 |   1 |

## Per-system runtime (total ms)

| System |   N |     Mean |   Median |      Std |      Min |      Max |
| :----- | --: | -------: | -------: | -------: | -------: | -------: |
| static |  90 |      169 |      150 |     41.2 |      124 |      309 |
| 2R     |  90 | 2.58e+04 | 2.03e+04 | 2.97e+04 |       38 | 1.24e+05 |
| 4R     |  90 | 4.18e+04 | 3.26e+04 | 3.42e+04 | 6.96e+03 | 1.35e+05 |

## Pairwise comparisons (paired, two-sided)

| A vs B       |   N | mean(B-A) |    t |      p_t |        W | p_w |  p_t_adj | p_w_adj | Cohen_dz |  r_rb | win_rate(B>A) | 95% CI mean(B-A) |
| :----------- | --: | --------: | ---: | -------: | -------: | --: | -------: | ------: | -------: | ----: | ------------: | ---------------- |
| static vs 2R |  90 |     0.252 | 9.76 | 9.84e-16 | 3.43e+03 |  NA | 1.97e-15 |      NA |     1.03 | 0.864 | [0.932,0.202] |
| static vs 4R |  90 |     0.269 | 11.3 | 6.67e-19 | 3.63e+03 |  NA |    2e-18 |      NA |     1.19 |     1 |    [ 1,0.222] |
| 2R vs 4R     |  90 |    0.0167 | 2.23 |   0.0282 |      440 |  NA |   0.0282 |      NA |    0.235 |     1 |   [ 1,0.0037] |

Notes:

- mean(B-A) > 0 means system B outperforms A on the chosen metric.
- Cohen_dz is effect size for paired samples (mean diff / sd diff).
- r_rb is the paired rank-biserial correlation; win_rate is fraction of bundles where B > A among non-ties.
- p-values adjusted with Holm-Bonferroni per test family (t-tests, Wilcoxon) across all pairs.

---

# Easy tier

# Bench Statistics (f1_edges)

## Per-system summary (over bundles)

| System |   N | Mean | Median | Std | Min | Max |
| :----- | --: | ---: | -----: | --: | --: | --: |
| static |  30 |    1 |      1 |   0 |   1 |   1 |
| 2R     |  30 |    1 |      1 |   0 |   1 |   1 |
| 4R     |  30 |    1 |      1 |   0 |   1 |   1 |

## Per-system runtime (total ms)

| System |   N |     Mean |   Median |      Std |      Min |      Max |
| :----- | --: | -------: | -------: | -------: | -------: | -------: |
| static |  30 |      132 |      128 |     12.2 |      124 |      192 |
| 2R     |  30 |     40.5 |       39 |     8.48 |       38 |       85 |
| 4R     |  30 | 9.68e+03 | 8.98e+03 | 3.39e+03 | 6.96e+03 | 2.58e+04 |

## Pairwise comparisons (paired, two-sided)

| A vs B       |   N | mean(B-A) |   t | p_t |   W | p_w | p_t_adj | p_w_adj | Cohen_dz | r_rb | win_rate(B>A) | 95% CI mean(B-A) |
| :----------- | --: | --------: | --: | --: | --: | --: | ------: | ------: | -------: | ---: | ------------: | ---------------- |
| static vs 2R |  30 |         0 |  NA |  NA |   0 |   1 |      NA |       1 |       NA |    0 |     [ 0.5, 0] |
| static vs 4R |  30 |         0 |  NA |  NA |   0 |   1 |      NA |       1 |       NA |    0 |     [ 0.5, 0] |
| 2R vs 4R     |  30 |         0 |  NA |  NA |   0 |   1 |      NA |       1 |       NA |    0 |     [ 0.5, 0] |

Notes:

- mean(B-A) > 0 means system B outperforms A on the chosen metric.
- Cohen_dz is effect size for paired samples (mean diff / sd diff).
- r_rb is the paired rank-biserial correlation; win_rate is fraction of bundles where B > A among non-ties.
- p-values adjusted with Holm-Bonferroni per test family (t-tests, Wilcoxon) across all pairs.

---

# Hard tier

# Bench Statistics (f1_edges)

## Per-system summary (over bundles)

| System |   N |  Mean | Median |    Std | Min |   Max |
| :----- | --: | ----: | -----: | -----: | --: | ----: |
| static |  30 | 0.539 |    0.5 | 0.0717 | 0.5 | 0.667 |
| 2R     |  30 |     1 |      1 |      0 |   1 |     1 |
| 4R     |  30 |     1 |      1 |      0 |   1 |     1 |

## Per-system runtime (total ms)

| System |   N |     Mean |  Median |      Std |      Min |      Max |
| :----- | --: | -------: | ------: | -------: | -------: | -------: |
| static |  30 |      160 |     148 |     27.9 |      141 |      242 |
| 2R     |  30 | 2.56e+04 | 2.2e+04 | 8.65e+03 | 1.53e+04 | 4.12e+04 |
| 4R     |  30 | 4.19e+04 | 3.5e+04 |  1.5e+04 | 2.58e+04 | 7.63e+04 |

## Pairwise comparisons (paired, two-sided)

| A vs B       |   N | mean(B-A) |    t |      p_t |   W | p_w |  p_t_adj | p_w_adj | Cohen_dz | r_rb | win_rate(B>A) | 95% CI mean(B-A) |
| :----------- | --: | --------: | ---: | -------: | --: | --: | -------: | ------: | -------: | ---: | ------------: | ---------------- |
| static vs 2R |  30 |     0.461 | 35.2 | 2.35e-25 | 465 |  NA | 4.71e-25 |      NA |     6.43 |    1 |    [ 1,0.433] |
| static vs 4R |  30 |     0.461 | 35.2 | 2.35e-25 | 465 |  NA | 2.35e-25 |      NA |     6.43 |    1 |    [ 1,0.433] |
| 2R vs 4R     |  30 |         0 |   NA |       NA |   0 |   1 |       NA |       1 |       NA |    0 |     [ 0.5, 0] |

Notes:

- mean(B-A) > 0 means system B outperforms A on the chosen metric.
- Cohen_dz is effect size for paired samples (mean diff / sd diff).
- r_rb is the paired rank-biserial correlation; win_rate is fraction of bundles where B > A among non-ties.
- p-values adjusted with Holm-Bonferroni per test family (t-tests, Wilcoxon) across all pairs.

---

# Veryhard tier

# Bench Statistics (f1_edges)

## Per-system summary (over bundles)

| System |   N |  Mean | Median |   Std |   Min |   Max |
| :----- | --: | ----: | -----: | ----: | ----: | ----: |
| static |  30 | 0.512 |    0.5 | 0.156 | 0.286 | 0.857 |
| 2R     |  30 | 0.809 |  0.762 | 0.165 |   0.5 |     1 |
| 4R     |  30 | 0.859 |  0.929 | 0.156 | 0.667 |     1 |

## Per-system runtime (total ms)

| System |   N |     Mean |   Median |      Std |      Min |      Max |
| :----- | --: | -------: | -------: | -------: | -------: | -------: |
| static |  30 |      214 |      209 |     25.4 |      195 |      309 |
| 2R     |  30 |  5.2e+04 | 3.44e+04 | 3.51e+04 | 1.88e+04 | 1.24e+05 |
| 4R     |  30 | 7.39e+04 | 7.19e+04 | 3.49e+04 | 2.99e+04 | 1.35e+05 |

## Pairwise comparisons (paired, two-sided)

| A vs B       |   N | mean(B-A) |    t |      p_t |   W | p_w |  p_t_adj | p_w_adj | Cohen_dz |  r_rb | win_rate(B>A) | 95% CI mean(B-A) |
| :----------- | --: | --------: | ---: | -------: | --: | --: | -------: | ------: | -------: | ----: | ------------: | ---------------- |
| static vs 2R |  30 |     0.296 | 6.29 | 7.21e-07 | 418 |  NA | 1.44e-06 |      NA |     1.15 | 0.724 | [0.862,0.201] |
| static vs 4R |  30 |     0.346 | 10.6 | 1.78e-11 | 465 |  NA | 5.35e-11 |      NA |     1.93 |     1 |    [ 1,0.283] |
| 2R vs 4R     |  30 |      0.05 | 2.34 |   0.0264 | 140 |  NA |   0.0264 |      NA |    0.427 |     1 |   [ 1,0.0111] |

Notes:

- mean(B-A) > 0 means system B outperforms A on the chosen metric.
- Cohen_dz is effect size for paired samples (mean diff / sd diff).
- r_rb is the paired rank-biserial correlation; win_rate is fraction of bundles where B > A among non-ties.
- p-values adjusted with Holm-Bonferroni per test family (t-tests, Wilcoxon) across all pairs.
