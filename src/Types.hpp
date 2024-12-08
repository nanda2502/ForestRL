#ifndef TYPES_HPP
#define TYPES_HPP

#include <vector>

using Repertoires = std::vector<std::vector<std::vector<size_t>>>; // repertoires[agent][tree][trait]

using Tree = std::vector<std::vector<size_t>>;

struct strategyExpectedValues {
    double payoff;
    double proximal;
};

struct Params {
    size_t num_agents = 1000;
    size_t num_trees = 100;
    size_t num_demonstrators = 10;
    size_t num_traits = 8;
    size_t num_iterations = 5e6;
    double lifetime_scale = 10.0;
    double learning_rate = 0.05;
    double temperature = 0.01;
    double innovation_rate = 0.05;
};

enum Strategy {
    Payoff,
    Proximal,
    Individual,
    None
};

struct TimeStepData {
    double propConstrained;
    size_t timestep;
    int strategy;
    int age;  
};

struct Result {
    std::vector<TimeStepData> timeSeriesData;
};

#endif // TYPES_HPP