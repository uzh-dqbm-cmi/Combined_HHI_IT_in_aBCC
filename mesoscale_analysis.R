library(readr)
library(tidyverse)
library(ggplot2)
library(dplyr)

vals <- read_csv("mesoscale_absolute_values.csv")
vals <- vals %>%
  dplyr::filter(
    !marker %in% c("IL-3") # always 0
  )
rename_dict <- c(
  "I-TAC"  = "CXCL11",
  "IP-10"   = "CXCL10",
  "MIG"    = "CXCL9")

vals$marker <- dplyr::recode(vals$marker, !!!rename_dict)

vals <- vals %>%
  dplyr::mutate(
    time = dplyr::recode(
      time,
      "V1"  = "w0",
      "V3"  = "w2",
      "V4"  = "w4",
      "V13" = "w26"
    )
  )
vals$time <- factor(vals$time, levels = c("w0", "w2", "w4", "w26"))

log2fc <- vals %>%
  dplyr::group_by(pat, marker) %>%
  dplyr::mutate(
    baseline = value[time == "w0"][1],
    log2fc = log2((value+0.0001) / (baseline+0.0001)) # adding to avoid 0-div
  ) %>%
  dplyr::filter(time != "w0") %>%
  dplyr::select(sample, plate, marker, pat, time, value, log2fc) %>%
  dplyr::ungroup()


pvals <- log2fc %>%
  dplyr::group_by(marker, time) %>%
  dplyr::summarise(
    n = sum(!is.na(log2fc)),
    p_value = if (n > 0) {
      wilcox.test(log2fc, mu = 0, alternative = "two.sided")$p.value
    } else {
      NA_real_
    },
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    p.adj = p.adjust(p_value, method = "BH")
  ) %>%
  dplyr::arrange(marker, time)

pvals <- pvals %>%
  dplyr::mutate(
    p.adj.signif = dplyr::case_when(
      is.na(p.adj)      ~ NA_character_,
      p.adj <= 0.0001   ~ "****",
      p.adj <= 0.001    ~ "***",
      p.adj <= 0.01     ~ "**",
      p.adj <= 0.05     ~ "*",
      TRUE              ~ "ns"
    )
  )

sig_df <- log2fc %>%
  group_by(marker) %>%
  summarise(y_pos = max(log2fc, na.rm = TRUE) * 1.12, .groups = "drop") %>%
  right_join(pvals, by = "marker")

ggplot(log2fc, aes(x = time, y = log2fc)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_boxplot(outlier.shape = NA, width = 0.7, alpha = 0.7) +
  geom_jitter(width = 0.12, size = 1.6, alpha = 0.7) +
  geom_text(
    data = sig_df,
    aes(x = time, y = y_pos, label = p.adj.signif),
    inherit.aes = FALSE,
    size = 5
  ) +
  facet_wrap(~marker, scales = "free_y") +
  labs(
    x = NULL,
    y = "log2 fold change vs w0",
  ) +
  theme_classic() +
  theme(
    strip.background = element_rect(fill = "grey95"),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )
marker_order <- c(
  "IL-12p70",
  "IL-18",
  "TNF-alpha",
  "CXCL9",
  "CXCL10",
  "CXCL11")

ggplot(log2fc[!(log2fc$marker %in% marker_order),], aes(x = time, y = log2fc)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_boxplot(outlier.shape = NA, width = 0.7, alpha = 0.7) +
  geom_jitter(width = 0.12, size = 1.6, alpha = 0.7) +
  geom_text(
    data = sig_df[!(sig_df$marker %in% marker_order),],
    aes(x = time, y = y_pos, label = p.adj.signif),
    inherit.aes = FALSE,
    size = 5
  ) +
  facet_wrap(~marker, scales = "free_y") +
  labs(
    x = NULL,
    y = "log2 fold change vs w0",
  ) +
  theme_classic() +
  theme(
    strip.background = element_rect(fill = "grey95"),
    axis.text.x = element_text(angle = 0, hjust = 0.5))

log2fc$marker <- factor(log2fc$marker, levels = marker_order)
pvals$marker   <- factor(pvals$marker, levels = marker_order)

log2fc6 <- log2fc %>%
  dplyr::filter(
    marker %in% c("IL-12p70", "IL-18", "TNF-alpha", "CXCL9", "CXCL10", "CXCL11")
  )

pvals6 <- log2fc6 %>%
  dplyr::group_by(marker, time) %>%
  dplyr::summarise(
    n = sum(!is.na(log2fc)),
    p_value = if (n > 0) {
      wilcox.test(log2fc, mu = 0, alternative = "two.sided")$p.value
    } else {
      NA_real_
    },
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    p.adj = p.adjust(p_value, method = "BH")
  ) %>%
  dplyr::arrange(marker, time)

pvals6 <- pvals6 %>%
  dplyr::mutate(
    p.adj.signif = dplyr::case_when(
      is.na(p.adj)      ~ NA_character_,
      p.adj <= 0.0001   ~ "****",
      p.adj <= 0.001    ~ "***",
      p.adj <= 0.01     ~ "**",
      p.adj <= 0.05     ~ "*",
      TRUE              ~ "ns"
    )
  )

sig_df6 <- log2fc6 %>%
  group_by(marker) %>%
  summarise(y_pos = max(log2fc, na.rm = TRUE) * 1.12, .groups = "drop") %>%
  right_join(pvals6, by = "marker")

sig_df6$marker <- factor(sig_df6$marker, levels = marker_order)
log2fc6$marker   <- factor(log2fc6$marker, levels = marker_order)

sig_df$y_pos_uniform <- 6.2
sig_df6$y_pos_uniform <- 6.2

ggplot(log2fc6, aes(x = time, y = log2fc)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_boxplot(outlier.shape = NA, width = 0.7, alpha = 0.7) +
  geom_jitter(width = 0.12, size = 1.6, alpha = 0.7) +
  geom_text(
    data = sig_df6,
    aes(x = time, y = y_pos_uniform, label = p.adj.signif),
    inherit.aes = FALSE,
    size = 5
  ) +
  facet_wrap(~marker, scales = "free_y") +
  labs(
    x = NULL,
    y = "log2FC from w0",
  ) +
  theme_classic() +
  theme(
    strip.background = element_rect(fill = "grey95"),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )


metadata <- read_table("mesoscale_pat_metadata.txt")
log2fc <- merge(log2fc, metadata, by.x = "pat", by.y = "ID", all.x = TRUE)
log2fc6 <- merge(log2fc6, metadata, by.x = "pat", by.y = "ID", all.x = TRUE)
log2fc6nh <- log2fc6 %>%
  dplyr::filter(Hematological_malignancy != "yes")

pvals6nh <- log2fc6nh %>%
  dplyr::group_by(marker, time) %>%
  dplyr::summarise(
    n = sum(!is.na(log2fc)),
    p_value = if (n > 0) {
      wilcox.test(log2fc, mu = 0, alternative = "two.sided")$p.value
    } else {
      NA_real_
    },
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    p.adj = p.adjust(p_value, method = "BH")
  ) %>%
  dplyr::arrange(marker, time)

pvals6nh <- pvals6nh %>%
  dplyr::mutate(
    p.adj.signif = dplyr::case_when(
      is.na(p.adj)      ~ NA_character_,
      p.adj <= 0.0001   ~ "****",
      p.adj <= 0.001    ~ "***",
      p.adj <= 0.01     ~ "**",
      p.adj <= 0.05     ~ "*",
      TRUE              ~ "ns"
    )
  )
sig_df6nh <- log2fc6nh %>%
  group_by(marker) %>%
  summarise(y_pos = max(log2fc, na.rm = TRUE) * 1.12, .groups = "drop") %>%
  right_join(pvals6nh, by = "marker")

sig_df6nh$marker <- factor(sig_df6nh$marker, levels = marker_order)
log2fc6nh$marker <- factor(log2fc6nh$marker, levels = marker_order)

sig_df6nh$y_pos_uniform <- 6.2

shape_map <- c(
  "CR" = 16,
  "PR" = 18,
  "SD" = 15,
  "PD" = 17,
  "Nan" = 1
)

pat_cols <- log2fc6 %>%
  dplyr::distinct(pat, patcols) %>%
  tibble::deframe()

ggplot(log2fc6nh, aes(x = time, y = log2fc)) +
  geom_hline(yintercept = 0, linetype = 2) +
  
  geom_boxplot(outlier.shape = NA, width = 0.7, alpha = 0.7) +
  
  geom_jitter(
    aes(color = patcols, shape = BoR26),
    width = 0.12,
    size = 2,
    alpha = 0.8
  ) +
  
  scale_color_identity(
    guide = "legend",
    breaks = unname(pat_cols),   # hex codes
    labels = names(pat_cols),    # patient IDs
    name = "Patient"
  ) +
  ylim(-1.1,6.23) + 
  scale_shape_manual(values = shape_map) +
  
  geom_text(
    data = sig_df6,
    aes(x = time, y = y_pos_uniform, label = p.adj.signif),
    inherit.aes = FALSE,
    size = 5
  ) +
  
  facet_wrap(~marker, scales = "free_y") +
  labs(
    x = NULL,
    y = "log2FC from w0",
    color = "Patient",
    shape = "Response"
  ) +
  
  theme_classic() +
  theme(
    strip.background = element_rect(fill = "grey95"),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )


ggplot(log2fc6, aes(x = time, y = log2fc)) +
  geom_hline(yintercept = 0, linetype = 2) +
  
  geom_boxplot(outlier.shape = NA, width = 0.7, alpha = 0.7) +
  
  geom_jitter(
    aes(color = patcols),
    width = 0.12,
    size = 1.6,
    alpha = 0.8
  ) +
  
  scale_color_identity(
    guide = "legend",
    breaks = unname(pat_cols),   # hex codes
    labels = names(pat_cols),    # patient IDs
    name = "Patient"
  ) +
  
  geom_text(
    data = sig_df6,
    aes(x = time, y = y_pos, label = p.adj.signif),
    inherit.aes = FALSE,
    size = 5
  ) +
  
  facet_wrap(~marker, scales = "free_y") +
  labs(
    x = NULL,
    y = "log2FC from w0",
    color = "Patient"
  ) +
  
  theme_classic() +
  theme(
    strip.background = element_rect(fill = "grey95"),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

ggplot(log2fc6nh, aes(x = time, y = log2fc)) +
  geom_hline(yintercept = 0, linetype = 2) +
  
  geom_boxplot(outlier.shape = NA, width = 0.7, alpha = 0.7) +
  
  geom_jitter(
    aes(color = patcols, shape = BoR26),
    width = 0.12,
    size = 2,
    alpha = 0.8
  ) +
  
  scale_color_identity(
    guide = "legend",
    breaks = unname(pat_cols),   # hex codes
    labels = names(pat_cols),    # patient IDs
    name = "Patient"
  ) +
  
  scale_shape_manual(values = shape_map) +
  
  geom_text(
    data = sig_df6,
    aes(x = time, y = y_pos, label = p.adj.signif),
    inherit.aes = FALSE,
    size = 5
  ) +
  
  facet_wrap(~marker, scales = "free_y") +
  labs(
    x = NULL,
    y = "log2FC from w0",
    color = "Patient",
    shape = "Response"
  ) +
  
  theme_classic() +
  theme(
    strip.background = element_rect(fill = "grey95"),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )
