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
    void printMeanStratEVs();
    void printMeanTreeEVs();
    void writeAveragesToCSV(const std::string& filename, const Params& params);
    std::vector<Strategy> chosenStrategy; // chosenStrategy[timestep]
    std::vector<treeType> chosenTree; // chosenTree[timestep]
    std::vector<size_t> chosenTreeIndex;
    std::vector<size_t> chosenAgent;
    std::vector<int> ageAtTimestep;
    std::vector<double> receivedPayoffs;
    int sumNone = 0;

private:
    Repertoires repertoires; // repertoires[agent][tree][trait]
    std::vector<StrategyExpectedValues> expectedValues; // expectedValues[agent] 
    std::vector<std::vector<double>> treeExpectedValues; // treeExpectedValues[tree]
    std::vector<size_t> lifetimes; // lifetimes[agent]
    std::vector<size_t> ages; // ages[agent]
    std::vector<std::vector<double>> payoffs; // payoffs[tree][trait]
    std::vector<size_t> uniqueIndices; // uniqueIndices[agent] 
    double payoffTotal;

    void initialize(const Params& params, const std::vector<Tree>& trees);
    std::vector<size_t> getUnknownTraits(size_t agentIndex, size_t treeIndex);
    std::vector<size_t> getUsefulDemonstrators(size_t focalAgent, size_t treeIndex, const std::vector<size_t>& demonstrators);
    std::vector<size_t> sampleDemoTraits(size_t focalAgent, size_t treeIndex, const std::vector<size_t>& demonstrators, std::mt19937& gen);
    size_t learnPayoffBased(size_t treeIndex, const std::vector<size_t>& usefulTraits, std::mt19937& gen);
    size_t learnProximal(size_t focalAgent, const std::vector<size_t>& usefulDemonstrators, size_t treeIndex, std::mt19937& gen);
    size_t learnProximalMarkov(size_t focalAgent, const std::vector<size_t>& usefulDemonstrators, size_t treeIndex, std::mt19937& gen);
    void update(size_t chosenTrait, size_t agentIndex, size_t treeIndex, Strategy strategy, const Tree& tree, const Params& params, std::mt19937& gen);
    bool isLearnable(size_t trait, size_t agentIndex, size_t treeIndex, const Tree& tree);
    size_t sampleUnexploredTree(size_t agentIndex, std::mt19937& gen);
    static bool isConstrained(size_t treeIndex);
    size_t sampleTreeFromEV(size_t agentIndex, const Params& params, std::mt19937& gen);
    std::vector<std::vector<double>> calculateAverages(const Params& params);
};

#endif // AGENTS_HPP