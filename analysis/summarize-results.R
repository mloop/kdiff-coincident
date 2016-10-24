# Purpose: calculate bias of estimator of D(h) under different simulation conditions

# Author: Matthew Shane Loop

# Preliminaries

## Packages
library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(boot)
library(cowplot)
library(devtools)

## User-defined functions
bias_function <- function(data, indices) {
  d <- data[indices, ]
  bias_b <- mean(d$d_hat) - max(d$true_d)
  return(bias_b)
}

percent_bias <- function(data) {
  percent_bias <- round(100 * (sum(data$d_hat - data$true_d) / sum(data$true_d)), digits = 0)
  return(percent_bias)
}

## Load data
results <- tbl_df(read.table(file = 'results.txt', sep = '\t', header = TRUE))
results <- arrange(results, condition, iteration)

## Bias


bias <- results %>%
  group_by(form, mismeasurement_probability, range) %>%
  nest() %>%
  mutate(
    bootstrap_replicates = map(.$data, ~boot(., statistic = bias_function, R = 1000)),
    bias_estimate = map_dbl(bootstrap_replicates, "t0"),
    bias_lower = map_dbl(bootstrap_replicates, ~boot.ci(., type = "norm")$normal[1, 2]),
    bias_upper = map_dbl(bootstrap_replicates, ~boot.ci(., type = "norm")$normal[1, 3]),
    percent_bias = map_dbl(.$data, ~percent_bias(.))
  ) %>%
  select(mismeasurement_probability, range, form, bias_estimate, bias_lower, bias_upper, percent_bias)

## Removed
bias_a <- bias %>%
  filter(form == 'true' | form == 'removed' | form == 'snapped') %>%
  mutate(form = factor(form, levels = c("true", "snapped", "removed"))) %>%
  mutate(
    form = ifelse(form == "removed", "Removed overlapping points",
                  ifelse(form == "snapped", "Overlapping points", "Original data set"))
  )

p_removed <- ggplot(data = filter(bias_a), aes(x = mismeasurement_probability, y = bias_estimate, color = form)) + 
  geom_pointrange(aes(ymin = bias_lower, ymax = bias_upper), position = position_dodge(width = 0.04)) + 
  geom_text(data = filter(bias_a, mismeasurement_probability == 0.5, form == "Overlapping points"), aes(x = mismeasurement_probability + 0.25, y = bias_estimate, label = paste(percent_bias, "% bias", sep = ""))) +
  scale_color_manual(values = c("red", "blue", "light green"), name = "") +
  facet_wrap(~ range, labeller = "label_both") +
  xlab("Proportion of overlapping points") +
  ylab("Bias") +
  theme_bw() +
  theme(text = element_text(size = 12)) +
  scale_x_continuous(limits = c(-0.05, 1), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5)) +
  ggtitle("Removing overlapping points prevents bias")

## Jittered
bias_b <- bias %>%
  filter(form == 'true' | form == 'jittered_0_05' | form == 'jittered_0_1' | form == 'jittered_0_15' | form == 'snapped') %>%
  mutate(
    form = factor(form, levels = c("true", "snapped", 'jittered_0_05', 'jittered_0_1', 'jittered_0_15')),
    form = ifelse(form == "jittered_0_05", "Jittered with radius 0.05",
                  ifelse(form == "jittered_0_1", "Jittered with radius 0.1", 
                         ifelse(form == "jittered_0_15", "Jittered with radius 0.15",
                                ifelse(form == "snapped", "Overlapping points", "Original data set"))))
  )

p_jittered <- ggplot(data = bias_b, aes(x = mismeasurement_probability, y = bias_estimate, color = form)) + 
  geom_pointrange(aes(ymin = bias_lower, ymax = bias_upper), position = position_dodge(width = 0.04)) + 
  #geom_line() + 
  geom_text(data = filter(bias_a, mismeasurement_probability == 0.5, form == "Overlapping points"), aes(x = mismeasurement_probability + 0.25, y = bias_estimate, label = paste(percent_bias, "% bias", sep = ""))) +
  scale_color_manual(values = c("orange", "purple", "light green", "red", "blue"), name = "") +
  facet_wrap(~ range, labeller = "label_both") +
  xlab("Proportion of overlapping points") +
  ylab("Bias") +
  theme_bw() +
  theme(text = element_text(size = 12)) +
  scale_x_continuous(limits = c(-0.05, 1), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5)) +
  ggtitle("Jittering overlapping points does not prevent bias")

p <- plot_grid(p_removed, p_jittered, ncol = 1, labels = c("A", "B"))

save_plot("bias-prevented.png",
       p,
       dpi = 300,
       base_width = 10,
       base_height = 10)

save_plot("bias-prevented.tiff",
          p,
          dpi = 300,
          base_width = 10,
          base_height = 10)

## Power

### Removed
power_removed <- results %>%
  group_by(form, mismeasurement_probability, range) %>%
  summarise(power = mean(reject)) %>%
  filter(form == 'true' | form == 'removed' | form == 'snapped') %>%
  select(mismeasurement_probability, form, range, power) %>%
  spread(form, power) %>%
  mutate(
    percent_power_loss_due_to_overlapping = round((true - snapped) / true * 100, digits = 1),
    percent_power_loss_prevented_by_removing = round((removed - snapped) / (true - snapped) * 100, digits = 0)
  ) %>%
  gather(form, power, removed, snapped, true) %>%
  mutate(
    form = ifelse(form == "removed", "Removed overlapping points",
                  ifelse(form == "snapped", "Overlapping points", "Original data set"))
  )

power_improvement <- results %>%
  group_by(form, mismeasurement_probability, range) %>%
  summarise(power = mean(reject)) %>%
  filter(form == 'true' | form == 'removed' | form == 'snapped') %>%
  spread(form, power)

p_removed <- ggplot(data = power_removed, aes(x = mismeasurement_probability, y = power * 100, color = form)) + 
  geom_pointrange(aes(ymin = (power - 1.96 * sqrt(power * (1 - power) / 1000)) * 100,
                      ymax = (power + 1.96 * sqrt(power * (1 - power) / 1000)) * 100)) + 
  geom_line() + 
  geom_text(data = filter(power_improvement, mismeasurement_probability == 0.5), aes(x = mismeasurement_probability + 0.05, y = removed * 100 - 6, label = paste("Removing points improved\n power by ", round(removed - snapped, digits = 2) * 100, " percentage points", sep = "")), color = "black", size = 3) +
  scale_color_manual(values = c("red", "blue", "light green"), name = "") +
  facet_wrap(~ range, labeller = "label_both") +
  xlab("Proportion of overlapping points") +
  ylab("Power (%)") +
  theme_bw() +
  theme(text = element_text(size = 12)) +
  scale_x_continuous(limits = c(-0.05, 1), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5)) +
  ggtitle("Removing overlapping points prevents power loss")

### Jittered
power_jittered <- results %>%
  group_by(form, mismeasurement_probability, range) %>%
  summarise(power = mean(reject)) %>%
  filter(form == 'true' | form == 'jittered_0_05' | form == 'jittered_0_1' | form == 'jittered_0_15' | form == 'snapped') %>%
  spread(form, power) %>%
  mutate(
    percent_power_loss_due_to_overlapping = round((true - snapped) / true * 100, digits = 1),
    percent_power_loss_prevented_by_jitter05 = round((jittered_0_05 - snapped) / (true - snapped) * 100, digits = 0),
    percent_power_loss_prevented_by_jitter1 = round((jittered_0_1 - snapped) / (true - snapped) * 100, digits = 0),
    percent_power_loss_prevented_by_jitter15 = round((jittered_0_15 - snapped) / (true - snapped) * 100, digits = 0)
  ) %>%
  gather(form, power, jittered_0_05, jittered_0_1, jittered_0_15, snapped, true) %>%
  mutate(
    form = ifelse(form == "jittered_0_05", "Jittered with radius 0.05",
                  ifelse(form == "jittered_0_1", "Jittered with radius 0.1", 
                         ifelse(form == "jittered_0_15", "Jittered with radius 0.15",
                                ifelse(form == "snapped", "Overlapping points", "Original data set"))))
  )

power_improvement <- results %>%
  group_by(form, mismeasurement_probability, range) %>%
  summarise(power = mean(reject)) %>%
  filter(form == 'true' | form == 'jittered_0_05' | form == 'jittered_0_1' | form == 'jittered_0_15' | form == 'snapped') %>%
  spread(form, power)

p_jittered <- ggplot(data = power_jittered, aes(x = mismeasurement_probability, y = round(power * 100, digits = 1), color = form)) + 
  geom_pointrange(aes(ymin = (power - 1.96 * sqrt(power * (1 - power) / 1000)) * 100,
                      ymax = (power + 1.96 * sqrt(power * (1 - power) / 1000)) * 100)) + 
  geom_line() + 
  geom_text(data = filter(power_improvement, mismeasurement_probability == 0.5), aes(x = mismeasurement_probability + 0.04, y = jittered_0_15 * 100 - 8, label = paste("Jittering with radius 0.15\n improved power by ", round(jittered_0_15 - snapped, digits = 2) * 100, " percentage points", sep = "")), color = "black", size = 3) +
  scale_color_manual(values = c("orange", "purple", "light green", "red", "blue"), name = "") +
  facet_wrap(~ range, labeller = "label_both") +
  xlab("Proportion of overlapping points") +
  ylab("Power (%)") +
  theme_bw() +
  theme(text = element_text(size = 12)) +
  scale_x_continuous(limits = c(-0.05, 1), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5)) +
  ggtitle("Jittering overlapping points did not prevent power loss")

p <- plot_grid(p_removed, p_jittered, ncol = 1, labels = c("A", "B"))

save_plot("power-losses-prevented.png",
          p,
          dpi = 300,
          base_width = 10,
          base_height = 10)
save_plot("power-losses-prevented.tiff",
          p,
          dpi = 300,
          base_width = 10,
          base_height = 10)

# Scatterplot
a <- results %>%
  group_by(mismeasurement_probability, range, form) %>%
  select(mismeasurement_probability, iteration, range, form, d_hat) %>%
  filter(mismeasurement_probability == 0.5, range == 0.15) %>%
  spread(form, d_hat)

snapped_vs_true <- ggplot(data = a, aes(x = true, y = snapped)) +
  geom_point(alpha = 0.5, color = "blue") +
  geom_point(aes(x = true, y = removed), alpha = 0.5, color = "red") +
  geom_smooth(method = "lm", color = "blue") +
  geom_smooth(aes(x = true, y = removed), method = "lm", color = "red") +
  theme_gray() +
  geom_abline(slope = 1)

removed_vs_true <- ggplot(data = a, aes(x = true, y = removed)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm") +
  theme_gray()

jittered_0_05_vs_true <- ggplot(data = a, aes(x = true, y = jittered_0_05)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm") +
  theme_gray()

jittered_0_1_vs_true <- ggplot(data = a, aes(x = true, y = jittered_0_1)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm") +
  theme_gray()

jittered_0_15_vs_true <- ggplot(data = a, aes(x = true, y = jittered_0_15)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm") +
  theme_gray()

p <- ggdraw() +
  draw_plot(snapped_vs_true, 0, 0.5, 0.5, 0.5) +
  draw_plot(removed_vs_true, 0.5, 0.5, 1, 1) +
  draw_plot(jittered_0_05_vs_true, 0, 0, 0.3, 0.5) +
  draw_plot(jittered_0_1_vs_true, 0.3, 0, 0.6, 0.5) +
  draw_plot(jittered_0_15_vs_true, 0.6, 0, 1, 0.5)

# Reproducibility
session_info()
