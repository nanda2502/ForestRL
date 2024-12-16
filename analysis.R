setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(dplyr)
library(ggplot2)
library(scales)
library(colorRamps)
library(ggbeeswarm)
rawdata <- read.csv("output.csv")


plotStrategyPref <- function(data) {
  plot <- ggplot(data, aes(x = TimeBin, y = PropPayoff, color = factor(PropConstrained))) +
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
  print(plot)
  return(plot)
}



rawdatafiltered <- rawdata %>%
  group_by(PropConstrained) %>%
  filter(Timestep > 0.9 * max(Timestep)) %>%
  ungroup() %>%
  mutate(AgeBin =  floor((Age - min(Age)) / (diff(range(Age)) / 20)) * (diff(range(Age)) / 20) + min(Age)) 
  
ggplot(rawdatafiltered, aes(x = AgeBin)) + 
  geom_bar() + 
  labs(y = "Number of Individuals", x = "Age Bin") +
  theme_minimal() 




plotAge <- function(data) {
  plot <- ggplot(data, aes(x = AgeBin, y = PropPayoff, color = factor(PropConstrained), group = factor(PropConstrained))) +
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
  print(plot)
  return(plot)
}

makeHeatmapFixed <- function(data, num_bins = 10) {
  data$PropConstrained_bins <- cut(data$PropConstrained, breaks = seq(0, 1, length.out = num_bins + 1), include.lowest = TRUE)
  data$PropPayoff_bins <- cut(data$PropPayoff, breaks = seq(0, 1, length.out = num_bins + 1), include.lowest = TRUE)
  
  heatmap_data <- as.data.frame(table(data$PropConstrained_bins, data$PropPayoff_bins))
  heatmap_data$Proportion <- heatmap_data$Freq / sum(heatmap_data$Freq)
  
  heatmap_data$PropConstrained_value <- sapply(as.character(heatmap_data$Var1), function(x) {
    range <- as.numeric(unlist(strsplit(gsub("\\(|\\]|\\[", "", x), ",")))
    mean(range)
  })
  
  heatmap_data$PropPayoff_value <- sapply(as.character(heatmap_data$Var2), function(x) {
    range <- as.numeric(unlist(strsplit(gsub("\\(|\\]|\\[", "", x), ",")))
    mean(range)
  })
  
  plot <- ggplot(heatmap_data, aes(x = PropConstrained_value, y = PropPayoff_value, fill = Proportion)) +
    geom_tile() +
    scale_fill_gradientn(colors = colorRampPalette(c("white", "blue", "green", "orange", "red"))(100)) +
    labs(x = "Proportion constrained", y = "Proportion payoff", fill = "Proportion") +
    theme_minimal()
  
  print(plot)
  return(plot)
}

makeHeatmapVariable <- function(data, num_bins = 10) {
  data$PropConstrainedChosen_bins <- cut(data$PropConstrainedChosen, breaks = seq(0, 1, length.out = num_bins + 1), include.lowest = TRUE)
  data$PropPayoff_bins <- cut(data$PropPayoff, breaks = seq(0, 1, length.out = num_bins + 1), include.lowest = TRUE)
  
  heatmap_data <- as.data.frame(table(data$PropConstrainedChosen_bins, data$PropPayoff_bins))
  heatmap_data$Proportion <- heatmap_data$Freq / sum(heatmap_data$Freq)
  
  heatmap_data$PropConstrainedChosen_value <- sapply(as.character(heatmap_data$Var1), function(x) {
    range <- as.numeric(unlist(strsplit(gsub("\\(|\\]|\\[", "", x), ",")))
    mean(range)
  })
  
  heatmap_data$PropPayoff_value <- sapply(as.character(heatmap_data$Var2), function(x) {
    range <- as.numeric(unlist(strsplit(gsub("\\(|\\]|\\[", "", x), ",")))
    mean(range)
  })
  
  plot <- ggplot(heatmap_data, aes(x = PropConstrainedChosen_value, y = PropPayoff_value, fill = Proportion)) +
    geom_tile() +
    scale_fill_gradientn(colors = colorRampPalette(c("white", "blue", "green", "orange", "red"))(100)) +
    labs(x = "Proportion constrained", y = "Proportion payoff", fill = "Proportion") +
    theme_minimal()
  
  print(plot)
  return(plot)
}



#### Fixed ####

cmd <- paste(
  "./runforest.exe",
  "--num_agents 1000",
  "--num_trees 100",
  "--num_demonstrators 10",
  "--num_traits 8",
  "--num_iterations 5000000",
  "--lifetime_scale 0.5",
  "--learning_rate 0.05",
  "--tree_learning_rate 0.00",
  "--temperature 0.01",
  "--tree_temperature 0.02",
  "--innovation_rate 0.05", 
  "--constrained_payoff_scale 1.0",
  "--start_prop 0.0",
  "--end_prop 1.0",
  "--prop_step 0.25",
  "--update_trees false"
)

system(cmd)

rawdata_fixed <- read.csv("output_fixed.csv.gz")


data1 <- rawdata_fixed %>%
  filter(Agent == 1)

##### Heatmap #####

data <- rawdata_fixed %>%
  group_by(PropConstrained) %>%
  filter(Timestep > 0.9 * max(Timestep)) %>%
  ungroup() %>%
  group_by(PropConstrained, Agent) %>%
  summarise(
    PropPayoff = mean(Strategy == 0),
    PropConstrainedChosen = mean(Tree == 1),
    .groups = 'drop'
  )

makeHeatmapFixed(data, 20)

##### Strategy Preference #####

data <- rawdata_fixed %>%
  mutate(TimeBin = floor(Timestep / 20000) * 20000) %>%
  group_by(PropConstrained, TimeBin) %>%
  summarise(
    PropPayoff = mean(Strategy == 0),
    PropConstrainedChosen = mean(Tree == 1),
    .groups = 'drop'
  )

plotStrategyPref(data)

##### Age ######

data <- rawdata_fixed %>%
  group_by(PropConstrained) %>%
  filter(Timestep > 0.9 * max(Timestep)) %>%
  ungroup() %>%
  mutate(AgeBin = cut(Age, quantile(Age, probs = seq(0, 1, length.out = 21)), include.lowest = TRUE, labels = FALSE)) %>%
  group_by(PropConstrained, AgeBin) %>%
  mutate(AgeBin = median(Age)) %>%
  ungroup() %>%
  group_by(PropConstrained, AgeBin) %>%
  summarise(
    PropPayoff = mean(Strategy == 0),
    .groups = 'drop')

plotAge(data)

##### Success & Payoffs ###### 
data <- rawdata_fixed %>%
  filter (
    Payoff > 0
  ) %>%
  group_by(PropConstrained) %>%
  filter(Timestep > 0.9 * max(Timestep)) %>%
  ungroup() %>%
  mutate(Success = ifelse(Payoff > 0, 1,0)) %>%
  #mutate(TimeBin = floor(Timestep / 100) * 100) %>%
  #group_by(PropConstrained, TimeBin) %>%
  group_by(PropConstrained, Agent, Tree, Strategy, Age) %>%
  summarise(
    PropPayoff = mean(Strategy == 0),
    PropSuccess = mean(Success == 1),
    AvgPayoff = mean(Payoff),
    .groups = 'drop'
  )

ggplot(data, aes(x = Age, y = PropSuccess, color = factor(Tree, labels = c("Flat", "Constrained")))) +
  geom_point(alpha = 0.6) + 
  geom_smooth(method = "loess", se = FALSE) + 
  labs(x = "Age", y = "Success Proportion", color = "Tree Type") +
  facet_wrap(~ factor(Strategy, labels = c("Payoff", "Proximal"))) +
  theme_minimal()

ggplot(data, aes(x = Age, y = AvgPayoff, color = factor(Tree, labels = c("Flat", "Constrained")))) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", se = FALSE) +
  labs(x = "Age", y = "Average Payoff", color = "Tree Type") +
  facet_wrap(~ factor(Strategy, labels = c("Payoff", "Proximal"))) +
  theme_minimal()

#### Variable ####
cmd <- paste(
  "./build/runforest.exe",
  "--num_agents 1000",
  "--num_trees 100",
  "--num_demonstrators 10",
  "--num_traits 8",
  "--num_iterations 5000000",
  "--lifetime_scale 0.4",
  "--learning_rate 0.1",
  "--tree_learning_rate 0.0",
  "--temperature 0.02",
  "--tree_temperature 0.02",
  "--innovation_rate 0.05", 
  "--constrained_payoff_scale 1.0",
  "--start_prop 0.5",
  "--end_prop 0.5",
  "--prop_step 0.0",
  "--updateTrees true"
)

system(cmd)

rawdata_variable <- read.csv("output_variable.csv")

##### Heatmap #####

data <- rawdata_variable %>%
  filter (PropConstrained == 0.5) %>%
  group_by(PropConstrained) %>%
  filter(Timestep > 0.9 * max(Timestep)) %>%
  ungroup() %>%
  #mutate(TimeBin = floor(Timestep / 100) * 100) %>%
  #group_by(PropConstrained, TimeBin) %>%
  group_by(PropConstrained, Agent) %>%
  summarise(
    PropPayoff = mean(Strategy == 0),
    PropConstrainedChosen = mean(Tree == 1),
    .groups = 'drop'
  )

makeHeatmapVariable(data, 50)

##### Strategy Preference #####

data <- rawdata_variable %>%
  mutate(TimeBin = floor(Timestep / 20000) * 20000) %>%
  group_by(PropConstrained, TimeBin) %>%
  summarise(
    PropPayoff = mean(Strategy == 0),
    PropConstrainedChosen = mean(Tree == 1),
    .groups = 'drop'
  )

plotStrategyPref(data)

##### Age ######


data <- rawdata_variable %>%
  group_by(PropConstrained) %>%
  filter(Timestep > 0.9 * max(Timestep)) %>%
  ungroup() %>%
  mutate(AgeBin = cut(Age, quantile(Age, probs = seq(0, 1, length.out = 11)), include.lowest = TRUE, labels = FALSE)) %>%
  group_by(PropConstrained, AgeBin) %>%
  mutate(AgeBin = median(Age)) %>%
  ungroup() %>%
  group_by(PropConstrained, AgeBin) %>%
  summarise(
    PropPayoff = mean(Strategy == 0),
    .groups = 'drop')

plotAge(data)

##### Success rate and payoffs ###### 
data <- rawdata_variable %>%
  filter (
    PropConstrained == 0.5
    ,Payoff > 0
    ) %>%
  group_by(PropConstrained) %>%
  filter(Timestep > 0.9 * max(Timestep)) %>%
  ungroup() %>%
  mutate(Success = ifelse(Payoff > 0, 1,0)) %>%
  #mutate(TimeBin = floor(Timestep / 100) * 100) %>%
  #group_by(PropConstrained, TimeBin) %>%
  group_by(PropConstrained, Agent, Tree, Strategy) %>%
  summarise(
    PropPayoff = mean(Strategy == 0),
    PropSuccess = mean(Success == 1),
    AvgPayoff = mean(Payoff),
    .groups = 'drop'
  )

ggplot(data, aes(x = factor(Tree, labels = c("Flat", "Constrained")), y = PropSuccess, color = factor(Tree, labels = c("Flat", "Constrained")))) +
  geom_beeswarm() +
  labs(x = "Tree Type", y = "Success Proportion", color = "Tree Type") +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
  facet_wrap(~ factor(Strategy, labels = c("Payoff", "Proximal")))

ggplot(data, aes(x = factor(Tree, labels = c("Flat", "Constrained")), y = AvgPayoff, color = factor(Tree, labels = c("Flat", "Constrained")))) +
  geom_beeswarm() +
  labs(x = "Tree Type", y = "Average Payoff", color = "Tree Type") +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
  facet_wrap(~ factor(Strategy, labels = c("Payoff", "Proximal")))

#### Repertoires ####

repertoires <- read.csv("repertoires.csv")
repertoires <- repertoires %>%
  group_by(Tree) %>%
  mutate(sum_average = sum(Average)) %>%
  ungroup() %>%
  mutate(Tree = factor(Tree, levels = unique(Tree[order(sum_average)])))

ggplot(repertoires, aes(x = Tree, y = factor(Trait), fill = Average)) +
  geom_tile() +
  scale_fill_gradientn(colors = colorRampPalette(c("white", "blue", "green", "orange", "red"))(100)) +
  labs(x = "    Constrained Trees | Unconstrained Trees", y = "Trait", fill = "Average") +
  theme_minimal() + 
  theme(axis.text.x = element_blank(),  
        axis.ticks.x = element_blank())

tree_sums <- repertoires %>%
  group_by(Tree) %>%
  summarize(sum_average = sum(Average))

first_50_avg <- tree_sums %>%
  slice(1:50) %>%
  summarize(first_50_avg = mean(sum_average))

last_50_avg <- tree_sums %>%
  slice(51:100) %>%
  summarize(last_50_avg = mean(sum_average))

data.frame(first_50_avg, last_50_avg)



#### Convenience Calculations ######


softmax_strategy <- function(payoff, proximal) {
  temperature <- 0.02
  scaledPayoff <- payoff/temperature
  scaledProximal <- proximal/temperature
  total <- exp(scaledPayoff) + exp(scaledProximal)
  prob_payoff <- exp(scaledPayoff)/total
  return (prob_payoff)
}


softmax_strategy(1, 1.05)

softmax_trees <- function(firstEV) {
  vec <- c(firstEV, rep(1, 99))
  temperature <- 0.02
  scaledVec <- vec / temperature
  total <- sum(exp(scaledVec))
  prob_first_element <- exp(scaledVec[1]) / total
  return(prob_first_element)
}


softmax_trees(1.05)


make_payoffs <- function(num_nodes) {
  spacing <- 4/num_nodes
  payoffs <- rep(0,num_nodes)
  for (i in 2:(num_nodes)) {
    payoffs[i] <- spacing * i
  }
  payoffs
}
  
mean(make_payoffs(20))  

1999 %% 1000
