#ifndef AGENTS_HPP
#define AGENTS_HPP

#include "Types.hpp"
#include <cstddef>
#include <random>

class Agents {

public:
    Agents(const Params& params, const std::vector<Tree>& trees);
    void learn(const Params& params, const std::vector<Tree>& trees);
    double computeProportion();

private:
    Repertoires repertoires; // repertoires[agent][tree][trait]
    std::vector<strategyExpectedValues> expectedValues; // expectedValues[agent] 
    std::vector<size_t> lifetimes; // lifetimes[agent]
    std::vector<std::vector<double>> payoffs; // payoffs[tree][trait]
    std::vector<int> chosenStrategy; // chosenStrategy[timestep]

    void initialize(const Params& params, const std::vector<Tree>& trees);
    std::vector<size_t> getUnknownTraits(size_t agentIndex, size_t treeIndex);
    std::vector<size_t> getUsefulDemonstrators(size_t focalAgent, size_t treeIndex, const std::vector<size_t>& demonstrators);
    std::vector<size_t> sampleTraits(size_t focalAgent, size_t treeIndex, const std::vector<size_t>& demonstrators);
    size_t learnPayoffBased(size_t treeIndex, const std::vector<size_t>& usefulTraits, std::mt19937& gen);
    size_t learnProximal(size_t focalAgent, const std::vector<size_t>& usefulDemonstrators, size_t treeIndex, std::mt19937& gen);
    void update(size_t chosenTrait, size_t agentIndex, size_t treeIndex, Strategy strategy, const Tree& tree, const Params& params, std::mt19937& gen);
    bool isLearnable(size_t trait, size_t agentIndex, size_t treeIndex, const Tree& tree);
};

#endif // AGENTS_HPP