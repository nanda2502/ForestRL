#include "Agents.hpp"
#include "Types.hpp"
#include "Utils.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <random>
#include <vector>


Agents::Agents(const Params& params) {
    initialize(params);
}

void Agents::initialize(const Params& params) {

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

Strategy sampleStrategy(const strategyExpectedValues& expectedValues, std::mt19937& gen) {
    //softmax the two values and then use them as probabilities to sample a strategy
    double total = exp(expectedValues.payoff) + exp(expectedValues.proximal);
    double payoff_prob = exp(expectedValues.payoff) / total;
    std::uniform_real_distribution<> dis(0, 1);
    
    if (dis(gen) < payoff_prob) {
        return Payoff;
    } 
    
    return Proximal;
}

std::vector<size_t> Agents::sampleDemonstrators(size_t focalAgent, const Params& params, std::mt19937& gen) {
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

    for (size_t demonstrator : demonstrators) {
        const std::vector<size_t>& demoTraits = repertoires[demonstrator][treeIndex];


        for (size_t traitIndex = 0; traitIndex < demoTraits.size(); traitIndex++) {
            if (demoTraits[traitIndex] == 1 && focalTraits[traitIndex] == 0) {
                usefulTraits.emplace_back(traitIndex);
                break;
            }
        }
    }    
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
    for (size_t traitIndex : usefulTraits) {
        size_t trait = usefulTraits[traitIndex];
        traitWeights[traitIndex] = std::pow(payoffs[treeIndex][trait], SENSITIVITY);
        weightSum += traitWeights[traitIndex];   
    }
    //normalize
    std::transform(traitWeights.begin(), traitWeights.end(), traitWeights.begin(), 
        [weightSum](double weight) {
            return weight / weightSum;
        });
    
    std::discrete_distribution<> distribution(traitWeights.begin(), traitWeights.end());

    size_t sampledTraitIndex = distribution(gen);
    return usefulTraits[sampledTraitIndex];
}

size_t Agents::learnProximal(size_t focalAgent, const std::vector<size_t>& usefulDemonstrators, size_t treeIndex, std::mt19937& gen) {
    std::vector<double> demonstratorWeights(usefulDemonstrators.size());
    std::vector<size_t> sampledTraits(usefulDemonstrators.size());
    double weightSum = 0.0;

    for (size_t demonstrator : usefulDemonstrators) {
         const std::vector<size_t>& demoTraits = repertoires[demonstrator][treeIndex];
        std::vector <size_t> newTraits;
        for (size_t traitIndex = 0; traitIndex < demoTraits.size(); traitIndex++) {
            if (demoTraits[traitIndex] == 1 && repertoires[focalAgent][treeIndex][traitIndex] == 0) {
                newTraits.emplace_back(traitIndex);
            }
        }

        if (!newTraits.empty()) {
            demonstratorWeights[demonstrator] = std::pow(2.0, 1.0 - newTraits.size());
            std::uniform_int_distribution<size_t> traitDist(0, newTraits.size() - 1);
            sampledTraits[demonstrator] = newTraits[traitDist(gen)];
        } else {
            demonstratorWeights[demonstrator] = 0.0;
            sampledTraits[demonstrator] = 0;
        }
        weightSum += demonstratorWeights[demonstrator];        
    }
    // normalize the weights
    std::transform(demonstratorWeights.begin(), demonstratorWeights.end(), demonstratorWeights.begin(),
        [weightSum](double weight) {
            return weight / weightSum;
        });
    // sample a trait from the weighted distribution
    std::discrete_distribution<> distribution(demonstratorWeights.begin(), demonstratorWeights.end());
    size_t sampledDemonstratorIndex = distribution(gen);
    return sampledTraits[sampledDemonstratorIndex];
}

bool Agents::isLearnable(size_t trait, size_t agentIndex, size_t treeIndex, const Tree& tree) {
    const std::vector<size_t>& focalTraits = repertoires[agentIndex][treeIndex];

    for (size_t parent = 0; parent < tree[trait].size(); ++parent) {
        if (tree[trait][parent] == 1) {  // Assuming adjacency matrix uses 1 for presence of an edge
            if (focalTraits[parent] == 0) {  // Check if the parent trait is not learned
                return false;
            }
        }
    }
    return true;
}

void Agents::update(size_t chosenTrait, size_t agentIndex, size_t treeIndex, Strategy strategy, const Tree& tree, const Params& params) {

    if (!chosenTrait == 0) {
        //check if learnable
        if (isLearnable(chosenTrait, agentIndex, treeIndex, tree)) {
            repertoires[agentIndex][treeIndex][chosenTrait] = 1;
            expectedValues[agentIndex].proximal =
                (1 - params.learning_rate) * expectedValues[agentIndex].proximal + params.learning_rate * 1.0; 

            expectedValues[agentIndex].payoff = 
                (1 - params.learning_rate) * expectedValues[agentIndex].payoff + params.learning_rate * 1.0; 

        } else {
            expectedValues[agentIndex].proximal =
                (1 - params.learning_rate) * expectedValues[agentIndex].proximal + params.learning_rate * -1.0; 

            expectedValues[agentIndex].payoff = 
                (1 - params.learning_rate) * expectedValues[agentIndex].payoff + params.learning_rate * -1.0; 
        }
    }

    lifetimes[agentIndex] -= 1;
    if (lifetimes[agentIndex] == 0) {
        //reset the agent's repertoire
        repertoires[agentIndex][treeIndex] = std::vector<size_t>(tree.size(), 0);
        lifetimes[agentIndex] = sampleLifetime(params, agentIndex);
    }
}

// This function samples a single individual for learning. 
// Whether learning is done socially or by innovation depends 
// on the innovation rate Parameter.
void Agents::learn(const Params& params, const std::vector<Tree>& trees) {
    std::random_device rd;
    std::mt19937 gen(rd());

    auto treeIndex = sampleIndex(params.num_trees, gen);
    auto agentIndex = sampleIndex(params.num_agents, gen);
    auto traits = getUnknownTraits(agentIndex, treeIndex);
    size_t chosenTrait;    
    std::uniform_real_distribution<> learn_dis(0, 1);
    Strategy strategy;

    if (learn_dis(gen) < params.innovation_rate) {
        // innovation
        strategy = Individual;
        chosenTrait = sampleIndex(traits.size(), gen);
    } else {
        // social learning

        strategy = sampleStrategy(expectedValues[agentIndex], gen); 

        auto demonstrators = sampleDemonstrators(agentIndex, params, gen);
        
        if (strategy == Payoff) {
            // payoff learning
            auto usefulTraits = sampleTraits(agentIndex, treeIndex, demonstrators);
            chosenTrait = learnPayoffBased(treeIndex, usefulTraits, gen);
        } else {
            // proximal learning
            auto usefulDemonstrators = getUsefulDemonstrators(agentIndex, treeIndex, demonstrators);
            chosenTrait = learnProximal(agentIndex, usefulDemonstrators, treeIndex, gen);
        }
        
    }
    
    update(chosenTrait, agentIndex, treeIndex, strategy, trees[treeIndex], params);
}