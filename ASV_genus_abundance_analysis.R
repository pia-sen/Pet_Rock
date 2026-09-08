#!/usr/bin/env Rscript
# ASV-within-genus abundance analysis. Base R only (no extra packages needed).
#
# Inputs (expected in the working directory):
#   asv_abundance_tbl_long.csv   columns: index, timepoint, microbe_asv, abundance
#   ASVs_taxonomy.csv            columns: microbe_asv, phylum, class, order, family, genus, species
#
# Outputs written to the working directory:
#   genus_kruskal_results.csv
#   genus_evenness_results.csv
#   dominant_asv_shuffle_results.csv
#   genus_temporal_results.csv
#   genus_temporal_pairwise.csv

# 1. Load and merge the two tables.
abund <- read.csv("asv_abundance_tbl_long.csv", row.names = 1, stringsAsFactors = FALSE)
tax   <- read.csv("ASVs_taxonomy.csv", fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)

# The two files use different ASV ID formats: "ASV1" vs "ASV_1". Normalize before merging.
abund$asv <- gsub("^ASV", "ASV_", abund$microbe_asv)
names(tax)[names(tax) == "microbe_asv"] <- "asv"

df <- merge(abund, tax[, c("asv", "genus")], by = "asv")

# Chronological order, not alphabetical: Nov -> Dec -> Jan -> Feb.
samples <- c("Pet_16S_1117_S249", "Pet_16S_1230_S250", "Pet_16S_0125_S251", "Pet_16S_0207_S252")
stopifnot(setequal(samples, unique(df$timepoint)))

# 2. Reshape to one row per ASV, one column per sample. Keep only ASVs with a
# genus assignment and complete data across all 4 samples.
wide <- reshape(df[, c("asv", "genus", "timepoint", "abundance")],
                 idvar = c("asv", "genus"), timevar = "timepoint", direction = "wide")
names(wide) <- sub("^abundance\\.", "", names(wide))
wide <- wide[!is.na(wide$genus), ]
wide <- wide[complete.cases(wide[, samples]), ]

cat(sprintf("Genus-classified ASVs with complete data: %d\n", nrow(wide)))
cat(sprintf("Genera represented: %d\n\n", length(unique(wide$genus))))

# 3. Per-genus significance test: do the ASVs in a genus differ in abundance
# across the 4 samples, or are they interchangeable? Kruskal-Wallis, plus
# epsilon-squared as an effect size since p-values alone don't show magnitude.
epsilon_squared <- function(H, n_asv, N) (H - n_asv + 1) / (N - n_asv)

genera <- sort(unique(wide$genus))
kw_rows <- list()

for (g in genera) {
  sub <- wide[wide$genus == g, samples]
  n_asv <- nrow(sub)

  if (n_asv < 2) {
    kw_rows[[g]] <- data.frame(genus = g, n_ASVs = n_asv, statistic = NA, p_value = NA, epsilon_sq = NA)
    next
  }

  values_per_asv <- as.list(as.data.frame(t(sub)))  # one vector of 4 values per ASV
  kw <- kruskal.test(values_per_asv)
  eps2 <- epsilon_squared(unname(kw$statistic), n_asv, n_asv * length(samples))

  kw_rows[[g]] <- data.frame(genus = g, n_ASVs = n_asv,
                              statistic = unname(kw$statistic), p_value = kw$p.value,
                              epsilon_sq = eps2)
}

kw_res <- do.call(rbind, kw_rows)
tested <- kw_res[!is.na(kw_res$p_value), ]

# Bonferroni is the conservative choice here: it controls the family-wise
# error rate across all 41 tested genera, rather than the expected proportion
# of false positives (FDR). FDR is kept for reference only.
tested$p_bonferroni <- p.adjust(tested$p_value, method = "bonferroni")
tested$p_fdr_bh <- p.adjust(tested$p_value, method = "BH")
tested$sig_bonferroni_0.05 <- tested$p_bonferroni < 0.05
tested <- tested[order(tested$p_value), ]

write.csv(tested, "genus_kruskal_results.csv", row.names = FALSE)

cat("Bonferroni-significant genera:\n")
print(tested[tested$sig_bonferroni_0.05, c("genus", "n_ASVs", "statistic", "p_value", "p_bonferroni", "epsilon_sq")])

# 4. Evenness per genus: is abundance spread across many ASVs, or
# concentrated in a few? Pielou's evenness loses resolution once a genus has
# hundreds of ASVs, so the Gini coefficient and "% of ASVs needed to reach
# 50% of total abundance" are reported as the more informative measures.
pielou_evenness <- function(x) {
  x <- x[x > 0]
  p <- x / sum(x)
  shannon_h <- -sum(p * log(p))
  shannon_h / log(length(x))
}

gini_coefficient <- function(x) {
  x <- sort(x)
  n <- length(x)
  cum <- cumsum(x)
  (n + 1 - 2 * sum(cum) / cum[n]) / n
}

evenness_rows <- list()

for (g in genera) {
  sub <- wide[wide$genus == g, samples]
  n_asv <- nrow(sub)
  if (n_asv < 3) next  # evenness isn't meaningful with fewer than 3 ASVs

  mean_abundance <- rowMeans(sub)
  sorted_desc <- sort(mean_abundance, decreasing = TRUE)
  cumulative_share <- cumsum(sorted_desc) / sum(sorted_desc)
  n_for_half <- which(cumulative_share >= 0.5)[1]

  evenness_rows[[g]] <- data.frame(genus = g, n_ASVs = n_asv,
                                    pielou_evenness = pielou_evenness(mean_abundance),
                                    gini = gini_coefficient(mean_abundance),
                                    n_ASVs_for_50pct_abundance = n_for_half,
                                    pct_ASVs_for_50pct_abundance = round(100 * n_for_half / n_asv, 1))
}

evenness_res <- do.call(rbind, evenness_rows)
evenness_res <- evenness_res[order(-evenness_res$pielou_evenness), ]
write.csv(evenness_res, "genus_evenness_results.csv", row.names = FALSE)

cat("\nMost even / diversified genera (n_ASVs >= 10, top 10):\n")
print(head(evenness_res[evenness_res$n_ASVs >= 10, ], 10))

# 5. Does the identity of the dominant ASV shuffle across samples, or is it
# consistently the same ASVs on top? Friedman test on the ASV x sample
# matrix within a genus, treating samples as blocks and ASVs as treatments.
# Null hypothesis: the 4 samples rank a genus's ASVs no more consistently
# than chance. Kendall's W = Friedman chi-squared / (n_samples * (n_ASVs - 1))
# is the effect size: W near 1 means the same ASVs dominate every sample;
# W near 0 means dominance is close to random from sample to sample.
#
# friedman.test() expects blocks (samples) as rows and treatments (ASVs) as
# columns, so the matrix is transposed relative to `wide`.
#
# Conservative note: friedman.test()'s p-value relies on a chi-squared
# approximation that is asymptotic in the number of BLOCKS. We only have
# 4 samples (blocks) -- large ASV counts do not fix this, since the
# approximation's validity depends on n_samples, not n_ASVs. So this
# function also computes an assumption-light permutation p-value directly
# from the data: under the null, an ASV's rank within a sample carries no
# information about its rank in another sample, so we repeatedly reshuffle
# each sample's ranks independently and rebuild the null distribution of
# the Friedman statistic from scratch. That permutation p-value, not the
# asymptotic one, is used for the significance calls below.
kendalls_w <- function(genus_name, min_asv = 3, n_perm = 5000, seed = 1) {
  sub <- wide[wide$genus == genus_name, samples]
  n_asv <- nrow(sub)
  if (n_asv < min_asv) return(NULL)

  sample_by_asv <- t(as.matrix(sub))
  n_samples <- nrow(sample_by_asv)
  n_asvs <- ncol(sample_by_asv)

  friedman_result <- friedman.test(sample_by_asv)
  chisq_obs <- unname(friedman_result$statistic)
  W <- chisq_obs / (n_samples * (n_asvs - 1))

  friedman_statistic_from_ranks <- function(rank_sums) {
    (12 / (n_samples * n_asvs * (n_asvs + 1))) * sum(rank_sums^2) - 3 * n_samples * (n_asvs + 1)
  }

  set.seed(seed)
  ranks <- t(apply(sample_by_asv, 1, rank))     # rank ASVs within each sample
  observed_stat <- friedman_statistic_from_ranks(colSums(ranks))
  permuted_stats <- replicate(n_perm, {
    permuted_ranks <- t(apply(ranks, 1, sample))  # reshuffle rank-to-ASV assignment, within each sample
    friedman_statistic_from_ranks(colSums(permuted_ranks))
  })
  p_permutation <- (1 + sum(permuted_stats >= observed_stat)) / (n_perm + 1)

  data.frame(genus = genus_name, n_ASVs = n_asvs, n_samples = n_samples,
             friedman_chisq = chisq_obs, df = unname(friedman_result$parameter),
             p_asymptotic = friedman_result$p.value, p_permutation = p_permutation,
             kendalls_W = W)
}

shuffle_all <- do.call(rbind, lapply(genera, kendalls_w))

# Bonferroni across every genus actually tested here, exactly as done for
# the genus-significance test above -- the permutation p-value is the
# conservative choice given only 4 samples, so it is what gets corrected
# and used for the significance call, not the asymptotic p-value.
shuffle_all$p_permutation_bonferroni <- p.adjust(shuffle_all$p_permutation, method = "bonferroni")
shuffle_all$reproducible_ranking_bonferroni_0.05 <- shuffle_all$p_permutation_bonferroni < 0.05
shuffle_all <- shuffle_all[order(shuffle_all$kendalls_W), ]
write.csv(shuffle_all, "dominant_asv_shuffle_results.csv", row.names = FALSE)

significant_genera <- tested$genus[tested$sig_bonferroni_0.05 & !is.na(tested$sig_bonferroni_0.05)]
shuffle_significant <- shuffle_all[shuffle_all$genus %in% significant_genera, ]

cat("\nRank-concordance (shuffle) results for the Bonferroni-significant genera:\n")
print(shuffle_significant)

cat("\nMost shuffled genera overall (lowest Kendall's W, n_ASVs >= 10):\n")
print(head(shuffle_all[shuffle_all$n_ASVs >= 10, ], 8))

cat("\nLeast shuffled genera overall (highest Kendall's W, n_ASVs >= 10):\n")
print(tail(shuffle_all[shuffle_all$n_ASVs >= 10, ], 8))

# 6. Does a genus's total abundance change across the 4 timepoints? Only run
# for genera already significant in step 3 -- if a genus's ASVs aren't even
# established to differ from each other, testing whether the genus's total
# changes over time adds another comparison without a solid basis for it.
# Friedman test again, but transposed relative to step 5: ASVs are blocks
# and samples are treatments, since the question here is whether abundance
# differs by sample, not whether ASV rank is reproducible.
temporal_change <- function(genus_name) {
  sub <- wide[wide$genus == genus_name, samples]  # ASVs (rows/blocks) x samples (cols/treatments)
  fried <- friedman.test(as.matrix(sub))
  data.frame(genus = genus_name, n_ASVs = nrow(sub),
             friedman_chisq = unname(fried$statistic), df = unname(fried$parameter),
             p_value = fried$p.value)
}

temporal_all <- do.call(rbind, lapply(significant_genera, temporal_change))
write.csv(temporal_all, "genus_temporal_results.csv", row.names = FALSE)

cat("\nTemporal abundance change (genera significant in step 3):\n")
print(temporal_all)

# Where that test is significant, localize which timepoints actually differ:
# paired Wilcoxon signed-rank tests (paired by ASV) across all 6 sample
# pairs, Bonferroni-corrected across those 6 comparisons per genus.
pairwise_temporal <- function(genus_name) {
  sub <- wide[wide$genus == genus_name, samples]
  pairs <- combn(samples, 2, simplify = FALSE)
  rows <- lapply(pairs, function(pr) {
    w <- wilcox.test(sub[[pr[1]]], sub[[pr[2]]], paired = TRUE)
    data.frame(genus = genus_name, sample_1 = pr[1], sample_2 = pr[2], p_value = w$p.value)
  })
  do.call(rbind, rows)
}

temporal_pairwise <- do.call(rbind, lapply(significant_genera, pairwise_temporal))
temporal_pairwise$p_bonferroni <- ave(temporal_pairwise$p_value, temporal_pairwise$genus,
                                       FUN = function(p) p.adjust(p, method = "bonferroni"))
write.csv(temporal_pairwise, "genus_temporal_pairwise.csv", row.names = FALSE)

cat("\nPairwise timepoint comparisons (Bonferroni-corrected across 6 pairs per genus):\n")
print(temporal_pairwise)

# Interpretation caveat (not just a statistical footnote): these are
# relative abundances that sum to ~1 within each sample, so ASVs are not
# statistically independent of one another -- one ASV's share rising
# mechanically pulls others down. The shuffle test therefore describes
# whether the RANK ORDER of a genus's ASVs is reproducible across these
# 4 specific samples, not an independent, causal claim about each ASV's
# biology, and with n = 4 samples it should be read as a pattern in this
# dataset rather than a generalizable result.
