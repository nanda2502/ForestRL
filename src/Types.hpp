#ifndef TYPES_HPP
#define TYPES_HPP

#include <vector>

//Access pattern is repertoires[agent][tree][trait]
using Repertoires = std::vector<std::vector<std::vector<size_t>>>;

using Tree = std::vector<std::vector<size_t>>;

struct strategyExpectedValues {
    double payoff;
    double proximal;
};

struct Params {
    size_t num_agents = 1000;
    size_t num_trees = 100;
    size_t num_demonstrators = 10;
    double learning_rate = 0.03;
    double innovation_rate = 0.01;
};

enum Strategy {
    Payoff,
    Proximal,
    Individual
};



#endif // TYPES_HPP