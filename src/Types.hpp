#ifndef TYPES_HPP
#define TYPES_HPP

#include <vector>

using Repertoires = std::vector<std::vector<std::vector<size_t>>>; // repertoires[agent][tree][trait]

using Tree = std::vector<std::vector<size_t>>;

struct strategyExpectedValues {
    double similarity;
    double proximal;
};

struct Params {
    size_t num_agents = 1000;
    size_t num_trees = 100;
    size_t num_demonstrators = 10;
    size_t num_traits = 10;
    size_t num_iterations = 4e6;
    double lifetime_scale = 0.5;
    double learning_rate = 0.03;
    double temperature = 0.5;
    double innovation_rate = 0.02;
};

enum Strategy {
    Similarity,
    Proximal,
    Individual,
    None
};

struct TimeStepData {
    double propConstrained;
    size_t timestep;
    int strategy;  
};

struct Result {
    std::vector<TimeStepData> timeSeriesData;
};

#endif // TYPES_HPP