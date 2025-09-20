This repository contains a minimal R workflow to double-check the metrics produced by the main Python benchmarking scripts.

Its purpose is to validate that reported values (e.g., Node/Edge P/R/F1, GCR, paired tests) are correct and reproducible, not artifacts of implementation.

Contents

- bench_results.jsonl – raw per-bundle outputs from the main pipeline.
- verify_bench.R – R script to recompute summaries, paired t-tests and Wilcoxon (Pratt), effect sizes, and bootstrap CIs.
- bench_report.md – verification report generated from the JSONL.

The report confirms or challenges the Python results, helping prevent misrepresentation due to coding mistakes.
