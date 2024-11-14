#ifndef AGENTS_HPP
#define AGENTS_HPP

#include "Types.hpp"
#include <cstddef>
#include <random>

class Agents {

public:
    Agents(const Params& params);
    void learn(const Params& params, const std::vector<Tree>& trees);

private:
    Repertoires repertoires; // repertoires[agent][tree][trait]
    std::vector<strategyExpectedValues> expectedValues; // expectedValues[agent] 
    std::vector<int> lifetimes; // lifetimes[agent]
    std::vector<std::vector<double>> payoffs; // payoffs[tree][trait]

    void initialize(const Params& params);
    std::vector<size_t> getUnknownTraits(size_t agentIndex, size_t treeIndex);
    static std::vector<size_t> sampleDemonstrators(size_t focalAgent, const Params& params, std::mt19937& gen);
    std::vector<size_t> getUsefulDemonstrators(size_t focalAgent, size_t treeIndex, const std::vector<size_t>& demonstrators);
    std::vector<size_t> sampleTraits(size_t focalAgent, size_t treeIndex, const std::vector<size_t>& demonstrators);
    size_t learnPayoffBased(size_t treeIndex, const std::vector<size_t>& usefulTraits, std::mt19937& gen);
    size_t learnProximal(size_t focalAgent, const std::vector<size_t>& usefulDemonstrators, size_t treeIndex, std::mt19937& gen);
    void update(size_t chosenTrait, size_t agentIndex, size_t treeIndex, Strategy strategy, const Tree& tree, const Params& params);
    bool isLearnable(size_t trait, size_t agentIndex, size_t treeIndex, const Tree& tree);
};

#endif // AGENTS_HPP