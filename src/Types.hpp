#ifndef TYPES_HPP
#define TYPES_HPP

#include <vector>

using Repertoires = std::vector<std::vector<std::vector<size_t>>>; // repertoires[agent][tree][trait]

using Tree = std::vector<std::vector<size_t>>;

struct StrategyExpectedValues {
    double payoff;
    double proximal;
};

struct Params {
    size_t num_agents = 1000;
    size_t num_trees = 100;
    size_t num_demonstrators = 10;
    size_t num_traits = 8;
    size_t num_iterations = 5e6;
    double lifetime_scale = 1.0;
    double learning_rate = 0.05;
    double tree_learning_rate = 0.02;
    double temperature = 0.02;
    double tree_temperature = 0.02;
    double innovation_rate = 0.05;
    double constrained_payoff_scale = 1.0;
    double start_prop = 0.0;
    double end_prop = 1.0;
    double prop_step = 0.25;
    bool update_trees = false;
};

enum Strategy {
    Payoff,     // Coded as 0
    Proximal,   // 1
    Individual, // 2 (not in final output file)
    None        // 3 (not in final output file)
};

enum treeType {
    Flat,
    Constrained
};

struct TimeStepData {
    double propConstrained;
    size_t timestep;
    int strategy;
    int age;
    int treeType; 
    size_t agent;
    size_t treeIndex;
    double payoff;
};

struct Result {
    std::vector<TimeStepData> timeSeriesData;
};

#endif // TYPES_HPP