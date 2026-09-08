#!/usr/bin/env Rscript
# ASV-within-genus abundance analysis. Base R only.
#
# Inputs (working directory):
#   asv_abundance_tbl_long.csv   columns: index, timepoint, microbe_asv, abundance
#   ASVs_taxonomy.csv            columns: microbe_asv, phylum, class, order, family, genus, species
#
# Outputs:
#   genus_kruskal_results.csv
#   genus_evenness_results.csv
#   dominant_asv_shuffle_results.csv

abund <- read.csv("asv_abundance_tbl_long.csv", row.names = 1, stringsAsFactors = FALSE)
tax   <- read.csv("ASVs_taxonomy.csv", fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)

# abund uses "ASV1", tax uses "ASV_1" -- normalize before merging
abund$asv <- gsub("^ASV", "ASV_", abund$microbe_asv)
names(tax)[names(tax) == "microbe_asv"] <- "asv"

df <- merge(abund, tax[, c("asv", "genus")], by = "asv")

samples <- c("Pet_16S_1117_S249", "Pet_16S_1230_S250", "Pet_16S_0125_S251", "Pet_16S_0207_S252")  # Nov-Dec-Jan-Feb
stopifnot(setequal(samples, unique(df$timepoint)))

wide <- reshape(df[, c("asv", "genus", "timepoint", "abundance")],
                 idvar = c("asv", "genus"), timevar = "timepoint", direction = "wide")
names(wide) <- sub("^abundance\\.", "", names(wide))
wide <- wide[!is.na(wide$genus), ]
wide <- wide[complete.cases(wide[, samples]), ]

cat(sprintf("Genus-classified ASVs with complete data: %d\n", nrow(wide)))
cat(sprintf("Genera represented: %d\n\n", length(unique(wide$genus))))

# per-genus Kruskal-Wallis: do a genus's ASVs differ in abundance across samples?
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

  values_per_asv <- as.list(as.data.frame(t(sub)))
  kw <- kruskal.test(values_per_asv)
  eps2 <- epsilon_squared(unname(kw$statistic), n_asv, n_asv * length(samples))

  kw_rows[[g]] <- data.frame(genus = g, n_ASVs = n_asv,
                              statistic = unname(kw$statistic), p_value = kw$p.value,
                              epsilon_sq = eps2)
}

kw_res <- do.call(rbind, kw_rows)
tested <- kw_res[!is.na(kw_res$p_value), ]

# Bonferroni across all tested genera; FDR kept for reference
tested$p_bonferroni <- p.adjust(tested$p_value, method = "bonferroni")
tested$p_fdr_bh <- p.adjust(tested$p_value, method = "BH")
tested$sig_bonferroni_0.05 <- tested$p_bonferroni < 0.05
tested <- tested[order(tested$p_value), ]

write.csv(tested, "genus_kruskal_results.csv", row.names = FALSE)

cat("Bonferroni-significant genera:\n")
print(tested[tested$sig_bonferroni_0.05, c("genus", "n_ASVs", "statistic", "p_value", "p_bonferroni", "epsilon_sq")])

# evenness per genus: spread across many ASVs vs concentrated in a few
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
  if (n_asv < 3) next

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

# does the same ASV dominate every sample, or does dominance shuffle?
# Friedman test on the ASV x sample matrix, samples as blocks, ASVs as treatments.
# Kendall's W = chisq / (n_samples * (n_ASVs - 1)): near 1 = same ASVs dominate
# every sample, near 0 = dominance is close to random across samples.
#
# friedman.test()'s asymptotic p-value is unreliable with only 4 blocks
# (samples), regardless of ASV count, so this also runs a permutation test:
# reshuffle each sample's ranks independently and rebuild the null
# distribution of the Friedman statistic. The permutation p-value, not the
# asymptotic one, drives the significance calls below.
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
  ranks <- t(apply(sample_by_asv, 1, rank))
  observed_stat <- friedman_statistic_from_ranks(colSums(ranks))
  permuted_stats <- replicate(n_perm, {
    permuted_ranks <- t(apply(ranks, 1, sample))
    friedman_statistic_from_ranks(colSums(permuted_ranks))
  })
  p_permutation <- (1 + sum(permuted_stats >= observed_stat)) / (n_perm + 1)

  data.frame(genus = genus_name, n_ASVs = n_asvs, n_samples = n_samples,
             friedman_chisq = chisq_obs, df = unname(friedman_result$parameter),
             p_asymptotic = friedman_result$p.value, p_permutation = p_permutation,
             kendalls_W = W)
}

shuffle_all <- do.call(rbind, lapply(genera, kendalls_w))

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

# Caveat: these are relative abundances summing to ~1 per sample, so ASVs
# aren't independent -- one ASV's share rising mechanically pulls others
# down. The shuffle test describes whether rank order is reproducible
# across these 4 specific samples, not an independent biological claim,
# and with n = 4 it's a pattern in this dataset, not a generalizable result.
