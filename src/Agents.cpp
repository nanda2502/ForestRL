#include "Agents.hpp"
#include "Types.hpp"
#include "Trees.hpp"
#include "Utils.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <random>
#include <vector>
#include <iostream>

extern int debug;

Agents::Agents(const Params& params, const std::vector<Tree>& trees) {
    repertoires = Repertoires(params.num_agents, std::vector<std::vector<size_t>>(params.num_trees, std::vector<size_t>(params.num_traits, 0)));
    expectedValues = std::vector<strategyExpectedValues>(params.num_agents, strategyExpectedValues{0.5, 0.5});
    lifetimes = std::vector<size_t>(params.num_agents);
    chosenStrategy = std::vector<int>();

    std::random_device rd;
    std::mt19937 gen(rd());

    auto payoffPair = generatePayoffs(params, gen);
    payoffs = payoffPair.first;
    payoffTotal = payoffPair.second;

    initialize(params, trees);
}

/* std::vector<double> makeConstrainedDist(const Params& params) {
    int inputnodes = params.num_traits - 2;
    double startingValue = 1.112 - 0.217 * inputnodes + 0.016 * inputnodes * inputnodes;
    constexpr double decay_rate = 0.5;
    std::vector<double> dist(params.num_traits -1);
    dist[0] = startingValue;
    for (size_t i = 1; i < params.num_traits - 1; ++i) {
        dist[i] = dist[i - 1] * decay_rate;
    }
    return dist;
} */

std::vector<double> makeConstrainedDist(const Params& params) {
    return std::vector<double>(params.num_traits - 1, 0.5);
}


std::vector<double> makeFlatDist(const Params& params) {
    return std::vector<double>(params.num_traits - 1, 1.0/(params.num_traits - 1));
}

size_t sampleLifetime(const Params& params, std::mt19937& gen) {
    double mean =  params.lifetime_scale * params.num_trees * params.num_traits;

    std::poisson_distribution<size_t> lifetime_dist(mean);

    return lifetime_dist(gen);
} 

void Agents::initialize(const Params& params, const std::vector<Tree>& trees) { 
    std::random_device rd;
    std::mt19937 gen(rd());

    for (size_t agent = 0; agent < params.num_agents; ++agent) {
        lifetimes[agent] = sampleLifetime(params, gen);
        for (size_t tree_idx = 0; tree_idx < params.num_trees; ++tree_idx) {
            repertoires[agent][tree_idx][0] = 1;

            auto tree = trees[tree_idx];
            std::vector<double> traitDist;
            if (tree[0][2] == 1) {
                traitDist = makeFlatDist(params);
                for (size_t trait = 1; trait < params.num_traits; ++trait) {
                    std::uniform_real_distribution<> dis(0, 1);
                    if (dis(gen) < traitDist[trait - 1]) {
                        repertoires[agent][tree_idx][trait] = 1;
                    } else break;
                }
            } else {
                traitDist = makeConstrainedDist(params);
                for (size_t trait = 1; trait < params.num_traits; ++trait) {
                    std::uniform_real_distribution<> dis(0, 1);
                    if (dis(gen) < traitDist[trait - 1]) {
                        repertoires[agent][tree_idx][trait] = 1;
                    } 
                }
            }
        }
    }
}

std::vector<size_t> Agents::getUnknownTraits(size_t agentIndex, size_t treeIndex) {
    std::vector<size_t> unknown_traits;
    const std::vector<size_t>& traits = repertoires[agentIndex][treeIndex];
    for (size_t i = 0; i < traits.size(); i++) {
        if (traits[i] == 0) {
            unknown_traits.push_back(i);
        }
    }
    return unknown_traits;
}    

Strategy sampleStrategy(const strategyExpectedValues& expectedValues, const Params& params, std::mt19937& gen) {
    if (debug >= 1) std::cout << "Sampling strategy" << '\n';
    if (debug >= 1) std::cout << "Payoff: " << expectedValues.payoff << '\n';
    if (debug >= 1) std::cout << "Proximal: " << expectedValues.proximal << '\n';
    
    // Adjusting for temperature in softmax
    double scaledPayoff = expectedValues.payoff / params.temperature;
    double scaledProximal = expectedValues.proximal / params.temperature;

    
    double total = exp(scaledPayoff) + exp(scaledProximal);
    double payoff_prob = exp(scaledPayoff) / total;
    std::uniform_real_distribution<> dis(0, 1);
    
    if (dis(gen) < payoff_prob) {
        return Payoff;
    } 
    
    return Proximal;
}

std::vector<size_t> sampleDemonstrators(size_t focalAgent, const Params& params, std::mt19937& gen) {
    std::vector<size_t> demonstrator_indices;
    //sample params.num_demonstrators agents, but not the focal agent
    for (size_t i = 0; i < params.num_demonstrators; i++) {
        auto index = sampleIndex(params.num_agents, gen);
        while (index == focalAgent) {
            index = sampleIndex(params.num_agents, gen);
        }
        demonstrator_indices.push_back(index);
    }
    return demonstrator_indices;
}

std::vector<size_t> Agents::sampleTraits(size_t focalAgent, size_t treeIndex, const std::vector<size_t>& demonstrators) {
    std::vector<size_t> usefulTraits;
    const std::vector<size_t>& focalTraits = repertoires[focalAgent][treeIndex];
    if (debug >= 1) std::cout << "Focal traits:\n" ;
    if (debug >= 1) printVector(focalTraits);

    for (size_t demonstrator : demonstrators) {
        const std::vector<size_t>& demoTraits = repertoires[demonstrator][treeIndex];


        for (size_t traitIndex = 0; traitIndex < demoTraits.size(); traitIndex++) {
            if (demoTraits[traitIndex] == 1 && focalTraits[traitIndex] == 0) {
                usefulTraits.emplace_back(traitIndex);
                break;
            }
        }
    }
    if (debug >= 1) std::cout << "Useful traits:\n" ;
    if (debug >= 1) printVector(usefulTraits);    
    return usefulTraits;
}

std::vector<size_t> Agents::getUsefulDemonstrators(size_t focalAgent, size_t treeIndex, const std::vector<size_t>& demonstrators) {
    std::vector<size_t> usefulDemonstrators;
    const std::vector<size_t>& focalTraits = repertoires[focalAgent][treeIndex];

    for (size_t demonstrator : demonstrators) {
        const std::vector<size_t>& demoTraits = repertoires[demonstrator][treeIndex];

        bool useful = false;
        for (size_t traitIndex = 0; traitIndex < demoTraits.size(); traitIndex++) {
            if (demoTraits[traitIndex] == 1 && focalTraits[traitIndex] == 0) {
                useful = true;
                break;
            }
        }
        if (useful) {
            usefulDemonstrators.emplace_back(demonstrator);
        }
    }
    return usefulDemonstrators;
}

size_t Agents::learnPayoffBased(size_t treeIndex, const std::vector<size_t>& usefulTraits, std::mt19937& gen) {
    constexpr double SENSITIVITY = 5.0;
    std::vector<double> traitWeights(usefulTraits.size(), 0.0);

    double weightSum = 0.0;
    for (size_t i = 0; i < usefulTraits.size(); i++) {
        size_t trait = usefulTraits[i];
        traitWeights[i] = std::pow(payoffs[treeIndex][trait], SENSITIVITY);
        weightSum += traitWeights[i];   
    }
    //normalize
    std::transform(traitWeights.begin(), traitWeights.end(), traitWeights.begin(), 
        [weightSum](double weight) {
            return weight / weightSum;
        });
    if (debug >= 1) std::cout << "Trait weights:\n";
    if (debug >= 1) printVector(traitWeights);

    std::discrete_distribution<> distribution(traitWeights.begin(), traitWeights.end());

    size_t sampledTraitIndex = distribution(gen);
    if (debug >= 1) std::cout << "Sampled trait: " << usefulTraits[sampledTraitIndex] << '\n';
    return usefulTraits[sampledTraitIndex];
}

size_t Agents::learnProximal(size_t focalAgent, const std::vector<size_t>& usefulDemonstrators, size_t treeIndex, std::mt19937& gen) {
    std::vector<double> demonstratorWeights(usefulDemonstrators.size());
    std::vector<size_t> sampledTraits(usefulDemonstrators.size());
    double weightSum = 0.0;
    if (debug >= 1) std::cout << "Finding new traits in demonstrators\n";
    for (size_t i = 0; i < usefulDemonstrators.size(); i++) {
        size_t demonstrator = usefulDemonstrators[i];
        const std::vector<size_t>& demoTraits = repertoires[demonstrator][treeIndex];
        std::vector<size_t> newTraits;
        
        for (size_t traitIndex = 0; traitIndex < demoTraits.size(); traitIndex++) {
            if (demoTraits[traitIndex] == 1 && repertoires[focalAgent][treeIndex][traitIndex] == 0) {
                newTraits.emplace_back(traitIndex);
            }
        }

        if (!newTraits.empty()) {
            if (debug >= 1) std::cout << "Number of new traits: " << newTraits.size() << '\n';
            demonstratorWeights[i] = std::pow(2.0, 1.0 - newTraits.size());
            std::uniform_int_distribution<size_t> traitDist(0, newTraits.size() - 1);
            sampledTraits[i] = newTraits[traitDist(gen)];
        } else {
            demonstratorWeights[i] = 0.0;
            sampledTraits[i] = 0;
        }
        weightSum += demonstratorWeights[i];        
    }

    if (weightSum == 0.0) {
        return 0;
    }

    // normalize the weights
    std::transform(demonstratorWeights.begin(), demonstratorWeights.end(), demonstratorWeights.begin(),
        [weightSum](double weight) {
            return weight / weightSum;
        });
    
    if (debug >= 1) std::cout << "Demonstrator weights:\n";
    if (debug >= 1) printVector(demonstratorWeights);
    if (debug >= 1) std::cout << "Sampled traits:\n"; 
    if (debug >= 1) printVector(sampledTraits);

    // sample a trait from the weighted distribution
    std::discrete_distribution<> distribution(demonstratorWeights.begin(), demonstratorWeights.end());
    size_t sampledDemonstratorIndex = distribution(gen);
    if (debug >= 1) std::cout << "Sampled trait" << sampledTraits[sampledDemonstratorIndex] << '\n';
    return sampledTraits[sampledDemonstratorIndex];
}

bool Agents::isLearnable(size_t trait, size_t agentIndex, size_t treeIndex, const Tree& tree) {
    const std::vector<size_t>& focalTraits = repertoires[agentIndex][treeIndex];

    for (size_t parent = 0; parent < tree[trait].size(); ++parent) {
        if (tree[parent][trait] == 1) {  
            if (focalTraits[parent] == 0) {  
                return false;
            }
        }
    }
    return true;
}

void Agents::update(size_t chosenTrait, size_t agentIndex, size_t treeIndex, Strategy strategy, const Tree& tree, const Params& params, std::mt19937& gen) {
    if(debug >= 1) std::cout << "Updating" << '\n';
    if (!(chosenTrait == 0)) {
        //check if learnable
        double feedback;
        if (isLearnable(chosenTrait, agentIndex, treeIndex, tree)) {
            repertoires[agentIndex][treeIndex][chosenTrait] = 1;
            feedback = payoffs[treeIndex][chosenTrait] / payoffTotal;
            //feedback = 1.0;
        } else {
            feedback = -1.0;
        }
        if (debug >= 1) std::cout << "Feedback: " << feedback << '\n';
        if (strategy == Payoff) {
            expectedValues[agentIndex].payoff =
                (1 - params.learning_rate) * expectedValues[agentIndex].payoff + params.learning_rate * feedback;
            chosenStrategy.push_back(Payoff);
        } else if (strategy == Proximal) {
            expectedValues[agentIndex].proximal =
                (1 - params.learning_rate) * expectedValues[agentIndex].proximal + params.learning_rate * feedback;
            chosenStrategy.push_back(Proximal);
        } 
    }

    lifetimes[agentIndex] -= 1;
    if (lifetimes[agentIndex] == 0) {
        if(debug >= 1) std::cout << "Agent reset" << '\n';
        //reset the agent's repertoire
        repertoires[agentIndex][treeIndex] = std::vector<size_t>(tree.size(), 0);
        repertoires[agentIndex][treeIndex][0] = 1;
        lifetimes[agentIndex] = sampleLifetime(params, gen);
    }
}

void Agents::learn(const Params& params, const std::vector<Tree>& trees) {
    std::random_device rd;
    std::mt19937 gen(rd());

    auto treeIndex = sampleIndex(params.num_trees, gen);
    auto agentIndex = sampleIndex(params.num_agents, gen);

    if (debug >= 1) {
        if (trees[treeIndex][0][2] == 1) {
            std::cout << "Flat tree" << '\n';
        } else {
            std::cout << "Constrained tree" << '\n';
        }   
    }
    
    size_t chosenTrait;    
    std::uniform_real_distribution<> learn_dis(0, 1);
    Strategy strategy;

    auto unknownTraits = getUnknownTraits(agentIndex, treeIndex);
    if(debug >= 1 && !unknownTraits.empty()) {
        std::cout << "Unknown traits:\n";
        printVector(unknownTraits);
    }
    if (unknownTraits.empty()) {
        if (debug >= 1) std::cout << "No unknown traits" << '\n';
        strategy = None;
        chosenTrait = 0;
    } else {
        if (learn_dis(gen) < params.innovation_rate) {
            // innovation
            if(debug >= 1) std::cout << "Innovation" << '\n';
            strategy = Individual;
            chosenTrait = unknownTraits[sampleIndex(unknownTraits.size(), gen)];
            if(debug >= 1) std::cout << "Innovation done" << '\n';
        } else {
            // social learning
            if(debug >= 1) std::cout << "Social learning" << '\n';
            strategy = sampleStrategy(expectedValues[agentIndex], params, gen); 
            auto demonstrators = sampleDemonstrators(agentIndex, params, gen);
            if (strategy == Payoff) {
                // payoff learning
                if(debug >= 1) std::cout << "Strategy: Payoff" << '\n';
                auto usefulTraits = sampleTraits(agentIndex, treeIndex, demonstrators);

                if (usefulTraits.empty()) {
                    strategy = None;
                    chosenTrait = 0;
                    sumNone++;
                } else {
                    chosenTrait = learnPayoffBased(treeIndex, usefulTraits, gen);
                }
            } else {
                // proximal learning
                if(debug >= 1) std::cout << "Strategy: Proximal" << '\n';
                auto usefulDemonstrators = getUsefulDemonstrators(agentIndex, treeIndex, demonstrators);

                if (usefulDemonstrators.empty()) {
                    strategy = None;
                    chosenTrait = 0;
                } else {
                    chosenTrait = learnProximal(agentIndex, usefulDemonstrators, treeIndex, gen);
                }
            }
            
        }
    }
    update(chosenTrait, agentIndex, treeIndex, strategy, trees[treeIndex], params, gen);
}

double Agents::computeProportion() {
    int numElements = chosenStrategy.size() * 0.05;
    auto start = chosenStrategy.end() - numElements;
    double sum = std::accumulate(start, chosenStrategy.end(), 0);
    return sum / numElements;
}