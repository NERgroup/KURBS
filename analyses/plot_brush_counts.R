################################################################################
# Plot processed brush sample counts
# Run after output/urchin_count_data/Step1_process_brush_data.R.
################################################################################

library(dplyr)
library(ggplot2)
library(lubridate)
library(readr)
library(scales)

processed_dir <- "/Volumes/enhydra/data/KURBS/recruitment_data/brushes/processed"
sample_type_data_file <- file.path(processed_dir, "sample_type_counts_processed.csv")
figure_dir <- "figures"

################################################################################
# Read plot-ready sample-type data

plot_dat <- readr::read_csv(sample_type_data_file, show_col_types = FALSE) %>%
  mutate(
    patch_type = factor(patch_type, levels = c("forest", "barren")),
    position = factor(position, levels = c("benthic", "suspended")),
    daily_counts_adj = daily_counts + 0.01
  ) %>%
  filter(
    !is.na(daily_counts),
    !is.na(patch_type),
    !is.na(position),
    is.finite(daily_counts)
  )

################################################################################
# Compare benthic and suspended counts within each patch type

urchin_rate_model <- glmmTMB::glmmTMB(
  total_urchins ~ patch_type * position + offset(log(soak_time)) +
    (1 | site) + (1 | collection_date) + (1 | site:collection_date),
  family = glmmTMB::nbinom2,
  data = plot_dat
)

position_emmeans <- emmeans::emmeans(
  urchin_rate_model,
  ~ position | patch_type
)

position_contrasts <- emmeans::contrast(
  position_emmeans,
  method = "pairwise",
  adjust = "none"
)

dodge_width <- 0.75

position_test_results <- as.data.frame(position_contrasts) %>%
  mutate(
    p_adj = p.adjust(p.value, method = "holm"),
    group1 = "benthic",
    group2 = "suspended",
    patch_x = match(as.character(patch_type), levels(plot_dat$patch_type)),
    xmin = patch_x - dodge_width / 4,
    xmax = patch_x + dodge_width / 4,
    p_label = if_else(
      p_adj < 0.001,
      "adjusted p < 0.001",
      paste0("adjusted p = ", formatC(p_adj, format = "f", digits = 2))
    )
  ) %>%
  left_join(
    plot_dat %>%
      group_by(patch_type) %>%
      summarise(
        # Bracket annotations use the log10 coordinates of the plotted y scale.
        y.position = log10(max(daily_counts_adj, na.rm = TRUE) * 1.35),
        .groups = "drop"
      ),
    by = "patch_type"
  )

################################################################################
# Figure styling

position_cols <- c(
  "benthic" = "#0072B2",
  "suspended" = "#D55E00"
)

patch_labs <- c(
  "forest" = "Forest",
  "barren" = "Barren"
)

position_labs <- c(
  "benthic" = "Benthic",
  "suspended" = "Suspended"
)

my_theme <- theme_classic(base_size = 13) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(color = "black", size = 13),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    strip.background = element_blank(),
    strip.text = element_text(size = 13, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 12, color = "gray30", hjust = 0),
    panel.grid = element_blank(),
    panel.border = element_blank()
  )

################################################################################
# Log-scale boxplot

p_box_log <- ggplot(
  plot_dat,
  aes(x = patch_type, y = daily_counts_adj, fill = position)
) +
  geom_boxplot(
    position = position_dodge(width = dodge_width),
    width = 0.62,
    outlier.shape = NA,
    alpha = 0.75,
    color = "black",
    linewidth = 0.45
  ) +
  geom_point(
    aes(color = position, group = position),
    position = position_jitterdodge(
      jitter.width = 0.12,
      dodge.width = dodge_width
    ),
    alpha = 0.75,
    size = 2.2,
    stroke = 0
  ) +
  scale_fill_manual(
    values = position_cols,
    labels = position_labs,
    name = "Position"
  ) +
  scale_color_manual(
    values = position_cols,
    labels = position_labs,
    name = "Position"
  ) +
  scale_x_discrete(labels = patch_labs) +
  scale_y_log10(
    labels = comma,
    expand = expansion(mult = c(0.04, 0.12))
  ) +
  labs(
    x = "Patch type",
    y = expression("Daily urchin recruits " ~ log[10] ~ "(count + 0.01)"),
    title = "Daily urchin recruitment by patch and position",
    subtitle = "Within-patch benthic vs suspended; Holm-adjusted p-values"
  ) +
  ggpubr::stat_pvalue_manual(
    position_test_results,
    label = "p_label",
    xmin = "xmin",
    xmax = "xmax",
    y.position = "y.position",
    tip.length = 0.01,
    bracket.size = 0.4,
    size = 3
  ) +
  guides(
    fill = guide_legend(override.aes = list(alpha = 0.75)),
    color = "none"
  ) +
  my_theme

p_box_log

################################################################################
# Histogram faceted by patch type

plot_dat_hist <- plot_dat %>%
  mutate(position = factor(position, levels = c("suspended", "benthic")))

p_hist <- ggplot(
  plot_dat_hist,
  aes(x = daily_counts_adj, fill = position)
) +
  geom_histogram(
    bins = 25,
    alpha = 0.85,
    position = "stack",
    color = "white",
    linewidth = 0.25
  ) +
  scale_fill_manual(
    values = position_cols,
    breaks = c("benthic", "suspended"),
    labels = position_labs,
    name = "Position"
  ) +
  scale_x_log10(
    labels = comma,
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
  facet_wrap(
    ~ patch_type,
    ncol = 1,
    labeller = labeller(patch_type = patch_labs)
  ) +
  labs(
    x = expression("Daily urchin recruits " ~ log[10] ~ "(count + 0.01)"),
    y = "Frequency",
    title = "Distribution of daily urchin recruitment"
  ) +
  my_theme

p_hist

################################################################################
# Weekly total recruits across all sampling

weekly_dat <- plot_dat %>%
  mutate(week = floor_date(collection_date, unit = "week", week_start = 1)) %>%
  group_by(week, position) %>%
  summarise(total_urchins = sum(total_urchins, na.rm = TRUE), .groups = "drop") %>%
  mutate(position = factor(position, levels = c("benthic", "suspended")))

p_weekly <- ggplot(
  weekly_dat,
  aes(x = week, y = total_urchins, color = position, group = position)
) +
  geom_line(linewidth = 1.1, alpha = 0.9) +
  geom_point(size = 2.6, alpha = 0.9) +
  scale_color_manual(
    values = position_cols,
    labels = position_labs,
    name = "Position"
  ) +
  scale_y_continuous(
    labels = comma,
    expand = expansion(mult = c(0, 0.12))
  ) +
  scale_x_date(
    date_breaks = "2 weeks",
    date_labels = "%b %d",
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  labs(
    x = "Week",
    y = "Total urchin recruits",
    title = "Weekly total urchin recruits across all sampling"
  ) +
  my_theme

p_weekly

################################################################################
# Save figures in the project

fs::dir_create(figure_dir, recurse = TRUE)

ggsave(
  file.path(figure_dir, "brush_boxplot_log.png"),
  plot = p_box_log,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(figure_dir, "brush_histogram_stacked.png"),
  plot = p_hist,
  width = 7,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(figure_dir, "brush_weekly_totals.png"),
  plot = p_weekly,
  width = 8,
  height = 5,
  dpi = 300
)
