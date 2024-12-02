setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(dplyr)
library(ggplot2)
library(scales)

rawdata <- read.csv("output.csv")
data <- rawdata %>%
  mutate(TimeBin = floor(Timestep / 20000) * 20000) %>%
  group_by(PropConstrained, TimeBin) %>%
  summarise(
    PropPayoff = mean(Strategy == 0),
    .groups = 'drop'
  )

ggplot(data, aes(x = TimeBin, y = PropPayoff, color = factor(PropConstrained))) +
  geom_line(linewidth = 1) +
  scale_color_viridis_d(name = "Prop. Constrained") +
  scale_x_continuous(
    name = "Time Step",
    labels = label_number(scale = 1e-6, suffix = "M"), # Convert to millions
    breaks = seq(0, max(data$TimeBin), by = 1000000)
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

rawdata <- read.csv("output.csv")

rawdatafiltered <- rawdata %>%
  group_by(PropConstrained) %>%
  filter(Timestep > 0.9 * max(Timestep)) %>%
  ungroup() %>%
  mutate(AgeBin =  floor((Age - min(Age)) / (diff(range(Age)) / 20)) * (diff(range(Age)) / 20) + min(Age)) 
  
ggplot(rawdatafiltered, aes(x = AgeBin)) + 
  geom_bar() + 
  labs(y = "Number of Individuals", x = "Age Bin") +
  theme_minimal() 


data <- rawdata %>%
  group_by(PropConstrained) %>%
  filter(Timestep > 0.9 * max(Timestep)) %>%
  ungroup() %>%
  mutate(AgeBin = cut(Age, quantile(Age, probs = seq(0, 1, length.out = 21)), include.lowest = TRUE, labels = FALSE)) %>%
  group_by(PropConstrained, AgeBin) %>%
  mutate(AgeBin = median(Age)) %>%
  ungroup() %>%
  group_by(PropConstrained, AgeBin) %>%
  summarise(PropPayoff = mean(Strategy == 0), .groups = 'drop')

ggplot(data, aes(x = AgeBin, y = PropPayoff, color = factor(PropConstrained))) +
  geom_point() +
  geom_line() +
  scale_color_viridis_d(name = "Prop. Constrained") +
  scale_x_continuous(
    name = "Age"
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

















softmax <- function(payoff, proximal) {
  temperature <- 0.02
  
  scaledPayoff <- payoff/temperature
  scaledProximal <- proximal/temperature
  total <- exp(scaledPayoff) + exp(scaledProximal)
  prob_payoff <- exp(scaledPayoff)/total
  return (prob_payoff)
}


softmax(1, 1.025)
