library(ggplot2)
library(dplyr)

data <- read.delim(
  "pat_recist_data.txt",
  header = TRUE,
  sep = "\t",
  na.strings = c("NA", "na", "")
)

metadata <- read.delim(
  "mesoscale_pat_metadata.txt",
  header = TRUE,
  sep = "\t",
  na.strings = c("NA", "na", "Nan", "NaN", "")
)

recist_order <- c("CR", "PR", "SD", "PD", "N/A")

data <- data %>%
  left_join(
    metadata %>% select(ID, patcols, patn, EoT),
    by = c("patient" = "ID")
  )

valid_data <- data %>%
  filter(
    !is.na(weeks),
    Staging %in% recist_order
  )

patient_info <- metadata %>%
  select(ID, patn, patcols, EoT) %>%
  distinct() %>%
  mutate(
    pat_number = as.numeric(gsub("[^0-9]", "", patn))
  )

patient_order <- patient_info %>%
  arrange(desc(pat_number)) %>%
  pull(patn)

best_resp_df <- data %>%
  filter(
    !is.na(weeks),
    weeks <= 26
  ) %>%
  mutate(
    # Death is considered N/A
    Staging_for_best_response = case_when(
      Staging == "Death" ~ "N/A",
      Staging %in% recist_order ~ Staging,
      TRUE ~ NA_character_
    ),
    Staging_for_best_response = factor(
      Staging_for_best_response,
      levels = recist_order,
      ordered = TRUE
    )
  ) %>%
  filter(!is.na(Staging_for_best_response)) %>%
  group_by(patn) %>%
  summarise(
    best_response_26 = min(Staging_for_best_response),
    .groups = "drop"
  )

last_assessment <- data %>%
  filter(!is.na(weeks)) %>%
  group_by(patn) %>%
  summarise(
    last_week = max(weeks),
    .groups = "drop"
  )

death_times <- data %>%
  filter(
    !is.na(weeks),
    Staging == "Death"
  ) %>%
  group_by(patn) %>%
  summarise(
    death_week = min(weeks),
    .groups = "drop"
  )

patient_bars <- patient_info %>%
  select(patn, EoT) %>%
  
  left_join(
    best_resp_df,
    by = "patn"
  ) %>%
  
  left_join(
    last_assessment,
    by = "patn"
  ) %>%

  left_join(
    death_times,
    by = "patn"
  ) %>%
  
  mutate(
    
    best_response_26 = case_when(
      !is.na(death_week) ~ "N/A",
      !is.na(best_response_26) ~ as.character(best_response_26),
      TRUE ~ "No response assessment"
    ),
    
    end_week = last_week,
    treatment_end = case_when(
      is.na(end_week) ~ NA_real_,
      is.na(EoT) ~ end_week,
      TRUE ~ pmin(EoT, end_week)
    ),
    
    followup_start = case_when(
      is.na(EoT) ~ end_week,
      TRUE ~ pmin(EoT, end_week)
    ),
    
    best_response_26 = factor(
      best_response_26,
      levels = c(
        "CR",
        "PR",
        "SD",
        "PD",
        "N/A"
      )
    )
  )

treatment_bars <- patient_bars %>%
  filter(
    !is.na(end_week),
    end_week > 0
  ) %>%
  mutate(
    x_start = 0,
    x_end = treatment_end
  ) %>%
  filter(x_end > x_start)


followup_bars <- patient_bars %>%
  filter(
    !is.na(end_week),
    !is.na(EoT),
    end_week > EoT
  ) %>%
  mutate(
    x_start = EoT,
    x_end = end_week
  ) %>%
  filter(x_end > x_start)

patient_bars <- patient_bars %>%
  mutate(
    patn = factor(
      patn,
      levels = patient_order
    )
  ) %>%
  arrange(best_response_26, patn) %>%
  mutate(
    patn = factor(
      patn,
      levels = unique(patn)
    )
  )

final_patient_order <- levels(patient_bars$patn)

treatment_bars <- treatment_bars %>%
  mutate(
    patn = factor(
      patn,
      levels = final_patient_order
    )
  )

followup_bars <- followup_bars %>%
  mutate(
    patn = factor(
      patn,
      levels = final_patient_order
    )
  )

recist_events <- valid_data %>%
  mutate(
    patn = factor(
      patn,
      levels = final_patient_order
    ),
    Staging_Method = paste(
      Staging,
      method,
      sep = "_"
    )
  )

death_events <- data %>%
  filter(
    !is.na(weeks),
    Staging == "Death"
  ) %>%
  mutate(
    patn = factor(
      patn,
      levels = final_patient_order
    )
  )

response_colors <- c(
  "CR" = "#99ff99",
  "PR" = "#e1ffe1",
  "SD" = "#ffffaa",
  "PD" = "#fcdcdb",
  "N/A" = "grey80")

darken_color <- function(col, factor = 0.90) {
  rgb_col <- col2rgb(col) / 255
  
  rgb(
    rgb_col[1] * factor,
    rgb_col[2] * factor,
    rgb_col[3] * factor
  )
}

lighten_color <- function(col, factor = 0.45) {
  
  rgb_col <- col2rgb(col) / 255
  
  rgb(
    rgb_col[1] + (1 - rgb_col[1]) * factor,
    rgb_col[2] + (1 - rgb_col[2]) * factor,
    rgb_col[3] + (1 - rgb_col[3]) * factor
  )
}

response_colors_dark <- sapply(
  response_colors,
  darken_color
)

response_colors_light <- sapply(
  response_colors,
  lighten_color
)

swimmer_plot <- ggplot() +

geom_rect(
  data = treatment_bars,
  aes(
    xmin = x_start,
    xmax = x_end,
    ymin = as.numeric(patn) - 0.325,
    ymax = as.numeric(patn) + 0.325,
    fill = best_response_26
  ),
  color = "grey60"
) +
geom_rect(
  data = followup_bars,
  aes(
    xmin = x_start,
    xmax = x_end,
    ymin = as.numeric(patn) - 0.325,
    ymax = as.numeric(patn) + 0.325,
    fill = best_response_26
  ),
  color = "grey60",
  alpha = 0.45
) +
geom_vline(
  xintercept = 26,
  color = "grey60",
  linetype = "dashed",
  alpha = 0.7,
  linewidth = 0.5
)+
geom_vline(
  xintercept = 32,
  color = "grey60",
  linetype = "dashed",
  alpha = 0.7,
  linewidth = 0.5
)+
geom_point(
  data = recist_events,
  aes(
    x = weeks,
    y = as.numeric(patn),
    shape = Staging_Method
  ),
  color = "black",
  size = 3.5,
  stroke = 1.1
) +
geom_point(
  data = death_events,
  aes(
    x = weeks,
    y = as.numeric(patn),
    shape = "Death"
  ),
  color = "black",
  size = 4,
  stroke = 1.2
) +
  scale_y_continuous(
    breaks = seq_along(final_patient_order),
    labels = final_patient_order,
    expand = c(0, 0)
  ) +
  
  scale_x_continuous(
    breaks = seq(0, 35, by = 5),
    expand = c(0, 0)
  ) +
  
  scale_fill_manual(
    values = response_colors_dark
  ) +
  
  scale_shape_manual(
    values = c(
      "CR_clinical" = 16,
      "PR_clinical" = 18,
      "SD_clinical" = 15,
      "PD_clinical" = 17,
      "CR_imaging"  = 1,
      "PR_imaging"  = 5,
      "SD_imaging"  = 0,
      "PD_imaging"  = 2,
      "Death"       = 4
    ),
    labels = c(
      "CR_clinical" = "CR (Clinical)",
      "CR_imaging" = "CR (Imaging)",
      "PR_clinical" = "PR (Clinical)",
      "PR_imaging" = "PR (Imaging)",
      "SD_clinical" = "SD (Clinical)",
      "SD_imaging" = "SD (Imaging)",
      "PD_clinical" = "PD (Clinical)",
      "PD_imaging" = "PD (Imaging)",
      "Death" = "Death"
    )
  ) +
  
  labs(
    x = "Weeks",
    y = "Patient",
    fill = "Best Response",
    shape = "RECIST Assessment & Event"
  ) +
  
  theme_minimal(base_size = 12) +
  
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(
      face = "bold",
      size = 9
    ),
    legend.position = "bottom",
    legend.box = "vertical",
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      color = "grey40"
    )
  )

print(swimmer_plot)

#### Scatter plot / spider plot ####
scatterdata <- data %>%
  filter(!is.na(change_BL))

scatterdata <- scatterdata %>%
  left_join(
    metadata %>% select(ID, patcols, patn),
    by = c("patient" = "ID")
  )
scatterdata <- scatterdata %>%
  mutate(
    patn = factor(patn, levels = patient_order)
  )

patient_colors <- setNames(
  patient_info$patcols,
  patient_info$patn
)

ggplot(
  scatterdata,
  aes(
    x = weeks,
    y = change_BL,
    color = patn,
    shape = Staging,
    group = patn
  )) +
  scale_shape_manual(
    values = c(
      "CR" = 16,  # filled circle
      "PR" = 18,  # filled triangle
      "SD" = 15,  # filled square
      "PD" = 17   # filled diamond
    )) +
#  ylim(-100,100) +
  geom_hline(
    yintercept = 0,
    color = "grey40",
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_hline(
    yintercept = 20,
    color = "red",
    linetype = "dashed",
    linewidth = 0.7
  ) +
  geom_hline(
    yintercept = -30,
    color = "darkgreen",
    linetype = "dashed",
    linewidth = 0.7
  ) +
  geom_vline(xintercept = 26, 
             color = "black", 
             linetype = "dashed", 
             alpha = 0.7, 
             linewidth = 0.5) +
  geom_line(
    linewidth = 0.8,
    alpha = 0.7
  ) +
  geom_point(
    size = 3,
    alpha = 0.9
  ) +
  scale_color_manual(
    values = patient_colors,
    breaks = patient_order,
    drop = FALSE
  ) +
  labs(
    x = "Weeks",
    y = "Change in target lesion from baseline (%)",
    color = "Patient"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )


time_to_response <- data %>%
  filter(
    !is.na(weeks),
    !is.na(Staging),
    Staging %in% c("PR", "CR")
  ) %>%
  group_by(patient) %>%
  summarise(
    time_to_response = min(weeks),
    first_response = Staging[which.min(weeks)],
    .groups = "drop"
  )

follow_up <- data %>%
  filter(!is.na(weeks)) %>%
  group_by(patient) %>%
  summarise(
    follow_up_time = max(weeks),
    .groups = "drop"
  )

median_follow_up <- median(
  follow_up$follow_up_time,
  na.rm = TRUE
)

median_time_to_response <- median(
  time_to_response$time_to_response,
  na.rm = TRUE
)

cat("Median time to response:", median_time_to_response, "weeks\n")
cat("Median follow-up:", median_follow_up, "weeks\n")