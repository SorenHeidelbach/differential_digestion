setwd("/user_data/sh/differential_digestion")
pacman::p_load(
  "data.table",
  "ggplot2",
  "magrittr",
  "ggside",
  "gt",
  "webshot2",
  "stringr",
  "scatterPlotMatrix",
  "GGally",
  "metR" # only for discretizing cotinous color scale
)


################################################################################
################################################################################
## All of this is related to the assembly provided by Mantas
## at "data/processed/assembly.fasta"
## PacBio + Illumina polished assembly
################################################################################
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
cov_files <- c(
  "data/processed/mapping/NP_MboI_R1041/NP_MboI_R1041.cov.tsv",
  "data/processed/mapping/NP_R1041_260bps_sub/NP_R1041_260bps_sub.cov.tsv",
  "data/processed/mapping/NP_R1041_400bps_sub/NP_R1041_400bps_sub.cov.tsv",
  "data/processed/mapping/NP_R1041_400bps_sub_7.7GB/NP_R1041_400bps_sub_7.7GB.cov.tsv"
)

# Load coverage profiles
cov <- lapply(
  cov_files,
  function(x) {
    sample_name <- basename(x) %>% gsub(".cov.tsv", "", .)
    dt <- fread(x)[, sample := sample_name]
    dt <- dt[, c("#rname", "endpos", "meandepth")]
    setnames(dt, c("meandepth", "#rname", "endpos"), c(sample_name, "contig", "length"))
    dt
  }
) %>%
  Reduce(function(x, y) merge(x, y, by = c("contig", "length"), all = TRUE), .)


# Plot digested vs. undigested
p_mboi_vs_400 <- ggplot(cov) +
  aes(x = NP_R1041_400bps_sub, y = NP_MboI_R1041, size = length) +
  geom_hline(yintercept = 1, linetype = 1, alpha = 0.35) +
  geom_vline(xintercept = 1, linetype = 1, alpha = 0.35) +
  geom_point(alpha = 0.1, shape = 16, col = "black", ) +
  scale_x_log10(limits = c(0.1, 1000)) +
  scale_y_log10(limits = c(0.1, 1000)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_size_continuous(range = c(1, 20), breaks = c(1e4, 1e5, 1e6, 1e7)) +
  theme_bw() +
  ylab("MboI digestion coverage") +
  xlab("Untreated coverage")
ggsave(
  filename = "results/figures/differential_coverage.png",
  plot = p_mboi_vs_400,
  width = 7,
  height = 6
)

# Plot undigested subset vs. undigested
p_mboi_vs_400 <- ggplot(cov) +
  aes(x = NP_R1041_400bps_sub, y = NP_R1041_400bps_sub_7.7GB, size = length) +
  geom_hline(yintercept = 1, linetype = 1, alpha = 0.35) +
  geom_vline(xintercept = 1, linetype = 1, alpha = 0.35) +
  geom_point(alpha = 0.1, shape = 16, col = "black", ) +
  scale_x_log10(limits = c(0.1, 1000)) +
  scale_y_log10(limits = c(0.1, 1000)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_size_continuous(range = c(1, 20), breaks = c(1e4, 1e5, 1e6, 1e7)) +
  theme_bw() +
  ylab("Subset coverage") +
  xlab("Full coverage")
ggsave(
  filename = "results/figures/differential_coverage_sub.png",
  plot = p_mboi_vs_400,
  width = 7,
  height = 6
)


# Plot digested vs. undigested subset
p_mboi_vs_400 <- ggplot(cov) +
  aes(x = NP_R1041_400bps_sub_7.7GB, y = NP_MboI_R1041, size = length) +
  geom_hline(yintercept = 1, linetype = 1, alpha = 0.35) +
  geom_vline(xintercept = 1, linetype = 1, alpha = 0.35) +
  geom_point(alpha = 0.1, shape = 16, col = "black", ) +
  scale_x_log10(limits = c(0.1, 1000)) +
  scale_y_log10(limits = c(0.1, 1000)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_size_continuous(range = c(1, 20), breaks = c(1e4, 1e5, 1e6, 1e7)) +
  theme_bw() +
  ylab("MboI digestion coverage") +
  xlab("Untreated coverage")
ggsave(
  filename = "results/figures/differential_coverage_sub.png",
  plot = p_mboi_vs_400,
  width = 7,
  height = 6
)

# Plot undigested (260bps) vs. undigested (400bps)
p_260_vs_400 <- ggplot(cov) +
  aes(x = NP_R1041_400bps_sub, y = NP_R1041_260bps_sub, size = length) +
  geom_hline(yintercept = 1, linetype = 1, alpha = 0.35) +
  geom_vline(xintercept = 1, linetype = 1, alpha = 0.35) +
  geom_point(alpha = 0.1, shape = 16, col = "black") +
  scale_x_log10(limits = c(0.1, 1000)) +
  scale_y_log10(limits = c(0.1, 1000)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_size_continuous(range = c(1, 20), breaks = c(1e4, 1e5, 1e6, 1e7)) +
  theme_bw() +
  ylab("Untreated coverage (260bps pore speed)") +
  xlab("Untreated coverage (400bps pore speed)")
ggsave(
  filename = "results/figures/differential_coverage_260_vs_400.png",
  plot = p_260_vs_400,
  width = 7,
  height = 6
)

# Plot coverage shift
cov[
    , diff_MboI_400 := NP_R1041_400bps_sub / NP_MboI_R1041
  ][
    , diff_MboI_400_weighted := diff_MboI_400 * (sum(NP_MboI_R1041) / sum(NP_R1041_400bps_sub))
  ][
    , diff_260_400 := NP_R1041_260bps_sub / NP_R1041_400bps_sub
  ]

cov_duplicated <- with(cov, cov[rep(1:length(length), as.integer(log(length, base = 2)))])
p_diff <- ggplot(cov) +
  aes(x = NP_R1041_400bps_sub, y = diff_MboI_400, size = length) +
  geom_point(alpha = 0.2, shape = 16, col = "black") +
  theme_bw() +
  geom_text(x = -1, y = -1.05, label = "Untreated : MboI", size = 3, vjust=0, hjust=0) +
  geom_text(x = -1, y = 0.05, label = "1 : 1", size = 3, vjust=0, hjust=0) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_text(x = -1, y = 1.05, label = "10 : 1", size = 3, vjust=0, hjust=0) +
  geom_hline(yintercept = 10, linetype = 2) +
  geom_text(x = -1, y = 2.05, label = "100 : 1", size = 3, vjust=0, hjust=0) +
  geom_hline(yintercept = 100, linetype = 2) +
  geom_text(x = -1, y = 3.05, label = "1000 : 1", size = 3, vjust=0, hjust=0) +
  geom_hline(yintercept = 1000, linetype = 2) +
  scale_x_log10(limits = c(0.1, NA)) +
  scale_y_log10(limits = c(0.1, NA)) +
  #stat_density2d(
  #  data = cov_duplicated, 
  #  aes(x = NP_R1041_400bps_sub, y = diff_MboI_400, alpha = ..level..),
  #  col = "red",
  #  n = 500,
  #  bins = 20
  #  ) +
  scale_size_continuous(range = c(1, 20), breaks = c(1e4, 1e5, 1e6, 1e7)) +
  scale_alpha_continuous(range = c(0.1, 1)) +
  ylab("Coverage ratio (Untreated / MboI)") +
  xlab("Coverage (Untreated)") +
  theme(ggside.panel.scale = 0.3)
ggsave(
  filename = "results/figures/coverage_shift.png",
  plot = p_diff,
  width = 7,
  height = 6
)





p_diff_density <- ggplot(cov) +
  aes(y = diff_MboI_400) +
  geom_density(aes(weight = length), fill = "gray") +
  scale_y_log10(limits = c(0.1, NA)) +
  theme_bw() +
    geom_text(x = 0, y = -1.05, label = "Untreated : MboI", size = 3, vjust=0, hjust=0) +
  geom_text(x = 0, y = 0.05, label = "1 : 1", size = 3, vjust=0, hjust=0) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_text(x = 0, y = 1.05, label = "10 : 1", size = 3, vjust=0, hjust=0) +
  geom_hline(yintercept = 10, linetype = 2) +
  geom_text(x = 0, y = 2.05, label = "100 : 1", size = 3, vjust=0, hjust=0) +
  geom_hline(yintercept = 100, linetype = 2) +
  geom_text(x = 0, y = 3.05, label = "1000 : 1", size = 3, vjust=0, hjust=0) +
  geom_hline(yintercept = 1000, linetype = 2) +
  labs(
    y = "Coverage ratio (Untreated / MboI)",
    x = "Density (weighted  with contig length)"
  )
ggsave(
  filename = "results/figures/diff_density.png",
  plot = p_diff_density,
  width = 7,
  height = 6
)



################################################################################
################################################################################
## All of this is related to the assembly generated with mmlong2
## at "mmlong2/tmp/polishing/asm_pol_lenfilt.fasta"
## generated from data/processed/NP/NP_R1041_400bps_sub.fastq
################################################################################
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

# General theme elements for all plots
plot_theme <- theme_bw(
  base_size = 12,
  base_family = "Helvetica"
)

################################################################################
# Data loading and prepping

# bins information from mmlong2
bins <- fread("mmlong2/results/mmlong2_bins.tsv")
bins[
  , quality := ifelse(Completeness >= 90 & Contamination <= 10, "HQ", "MQ")
  ][
    , Type := "mmlong2"
  ][
    , GC_Content := GC_Content * 100
  ]

# Load coverage profiles
coverage <- fread("mmlong2/tmp/binning/metabat_cov.tsv") %>%
  melt(id.vars = c("contigName", "contigLen", "totalAvgDepth")) %>%
  setnames(c("contig", "length", "mean_depth_all", "sample_and_type", "value"))
coverage[
    , ID := str_extract(sample_and_type, "^[^_]*")
  ][
    , value_type := fifelse(str_detect(sample_and_type, "var"), "var_depth", "mean_depth")
  ]
coverage <- coverage[
    , !"sample_and_type"
  ]

# Load sample IDs
coverage_sample_id <- fread("mmlong2/tmp/binning/reads.csv", header = FALSE) %>%
  setnames(c("sample_type", "path"))
coverage_sample_id[
    , sample := str_extract(path, "[^/]*$")
  ][
    , sample := str_replace(sample, ".fastq", "")
  ][
    , ID := as.character(seq_len(.N))
  ]
coverage_sample_id <- coverage_sample_id[
    , !"path"
  ]

# Merge coverage profiles with sample IDs
coverage <- coverage[coverage_sample_id, on = "ID"]
coverage <- dcast(
  coverage,
  contig + length + mean_depth_all ~ sample + value_type,
  value.var = "value"
)

# Load contig to bin mapping
contig_to_bin <- fread("mmlong2/tmp/binning/contig_bin.tsv", header = FALSE) %>%
  setnames(c("contig", "bin"))

# Add bin info to contigs
contigs <- merge(
    coverage,
    contig_to_bin,
    by = "contig",
    all.x = TRUE
  ) %>%
  merge(
    bins,
    by = "bin",
    all.x = TRUE
  )

# Load MboI coverage profile


################################################################################
# Plotting

# MAG quality plot (completeness vs contamination)
p_mag_qual <- ggplot(bins) +
  aes(x = Completeness, y = Contamination, size = Genome_Size) +
  geom_vline(xintercept = 90) +
  geom_point(alpha = 0.5) +
  plot_theme +
  scale_size_continuous(range = c(1, 10)) +
  geom_text(x = 91, y = 7.5, label = "HQ-MAGs", size = 5, hjust = 0) +
  geom_text(x = 89, y = 7.5, label = "MQ-MAGs", size = 5, hjust = 1)

ggsave(
  "results/figures/mmlong2_bins_quality.png",
  width = 8, height = 7,
  p_mag_qual
)

# MAG quality table
table_mag_quality <- bins[
  , .(.N), by = c("quality", "Type")
  ] %>%
  dcast(Type ~ quality) %>%
  `[`(
    , Total := HQ + MQ
  ) %>%
  gt() %>%
  tab_header(
    title = "Number of bins generated"
  ) %>%
  tab_style(
    style = list(
      cell_fill(color = "#b9ffb7")
      ),
    locations = cells_body(
      columns = Total,
      rows = Total == max(Total)
    )
  ) %>%
    tab_style(
    style = list(
      cell_fill(color = "#b9ffb7")
      ),
    locations = cells_body(
      columns = MQ,
      rows = MQ == max(MQ)
    )
  ) %>%
    tab_style(
    style = list(
      cell_fill(color = "#b9ffb7")
      ),
    locations = cells_body(
      columns = HQ,
      rows = HQ == max(HQ)
    )
  )
gtsave(
  table_mag_quality,
  "results/figures/mag_quality_table.png",
  width = 8, height = 7
)


# Coverage plots
plot_cov_vs_var <- ggplot(contigs) +
  aes(x = NP_R1041_400bps_sub_mean_depth, y = NP_R104_mean_depth) +
  geom_point(alpha = 0.2) +
  plot_theme +
  scale_x_log10() +
  scale_y_log10()
ggsave(
  "results/figures/contigs_cov_vs_var.png",
  plot_cov_vs_var,
  width = 8, height = 7
)




coverage_cols <- names(contigs)[grepl("mean_depth", names(contigs))]
plot_cols <- coverage_cols[c(2, 6, 10, 12, 14)]
plot_data <- cbind(
  contigs[
      , ..plot_cols
    ][
      , lapply(.SD, log10)
    ] %>%
  setnames(paste0(names(.), "_log10")),
  contigs[, c("length")]
) 
plot_data[plot_data == -Inf] <- NA

plot_data_no_na <- plot_data[complete.cases(plot_data),]

# Matrix plot
upper_density <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_density2d_filled(mapping = aes(fill = after_stat(level))) +
    scale_fill_discretised(low = "white", high = "#000000")
}
lower_scatter <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_smooth(method = "lm", col = "red") +
    geom_point(mapping = aes(size = length), alpha = 0.1, shape = 16) +
    scale_size_continuous(trans = "log10", range = c(0.1, 1))
}
diag_density <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_density(aes(weight = length), fill = "gray")
}
plot_coverage_profiles <- ggpairs(
  plot_data_no_na,
  columns = names(plot_data)[!grepl("length", names(plot_data))],
  columnLabels = c(
    "IL-201104_mean_depth" = "IL-201104 (log10)",
    "IL-201502_mean_depth" = "IL-201502 (log10)",
    "IL-201804_mean_depth" = "IL-201804 (log10)",
    "NP_R103_mean_depth" = "NP-R103 (log10)",
    "NP_R1041_400bps_sub_mean_depth" = "NP-R1041 (log10)"
  ),
  lower = list(continuous = lower_scatter),
  diag = list(continuous = diag_density)
) +
theme_bw()
ggsave(
  "results/figures/contigs_coverage_profiles_correlation.png",
  plot_coverage_profiles,
  width = 8, height = 7
)

plot_data[plot_data == NA] <- 0
plot_coverage_profiles <- ggpairs(
  plot_data,
  columns = names(plot_data)[!grepl("length", names(plot_data))],
  columnLabels = c(
    "IL-201104_mean_depth" = "IL-201104 (log10)",
    "IL-201502_mean_depth" = "IL-201502 (log10)",
    "IL-201804_mean_depth" = "IL-201804 (log10)",
    "NP_R103_mean_depth" = "NP-R103 (log10)",
    "NP_R1041_400bps_sub_mean_depth" = "NP-R1041 (log10)"
  ),
  lower = list(continuous = lower_scatter),
  diag = list(continuous = diag_density),
  upper = list(continuous = upper_density)
) +
theme_bw()
ggsave(
  "results/figures/contigs_coverage_profiles.png",
  plot_coverage_profiles,
  width = 8, height = 7
)
