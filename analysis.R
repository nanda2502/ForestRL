setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(dplyr)
library(ggplot2)
library(scales)


data <- read.csv("output.csv")

plot(data$PropConstrained, data$PropPayoff, xlab="Proportion of Constrained Trees", ylab="Proportion Payoff-based", type = "l", ylim = c(0,1))


data <- read.csv("output.csv") %>%
  mutate(TimeBin = floor(Timestep / 20000) * 20000) %>%
  group_by(PropConstrained, TimeBin) %>%
  summarise(
    PropPayoff = mean(Strategy == 1),
    .groups = 'drop'
  )

ggplot(data, aes(x = TimeBin, y = PropPayoff, color = factor(PropConstrained))) +
  geom_line(linewidth = 1) +
  scale_color_viridis_d(name = "Prop. Constrained") +
  scale_x_continuous(
    name = "Time Step",
    labels = label_number(scale = 1e-6, suffix = "M"), # Convert to millions
    breaks = seq(0, max(data$TimeBin), by = 500000)
  ) +
  scale_y_continuous(
    name = "Prop. Payoff",
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2)
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    text = element_text(size = 12),
    axis.text = element_text(size = 10)
  )


