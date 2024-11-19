#include "Agents.hpp"
#include "Types.hpp"
#include "Trees.hpp"
#include "Utils.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <random>
#include <vector>
#include <iostream>

extern int debug;

size_t sampleLifetime(const Params& params, std::mt19937& gen) {
    double mean = params.lifetime_scale * params.num_trees;

    std::poisson_distribution<size_t> lifetime_dist(mean);

    return lifetime_dist(gen);
} 

Agents::Agents(const Params& params, const std::vector<Tree>& trees) {
    repertoires = Repertoires(params.num_agents, std::vector<std::vector<size_t>>(params.num_trees, std::vector<size_t>(params.num_traits, 0)));
    expectedValues = std::vector<strategyExpectedValues>(params.num_agents, strategyExpectedValues{0.5, 0.5});
    lifetimes = std::vector<size_t>(params.num_agents);
    chosenStrategy = std::vector<int>();

    std::random_device rd;
    std::mt19937 gen(rd());

    auto payoffPair = generatePayoffs(params, gen);

    payoffs = payoffPair.first;
    totalPayoff = payoffPair.second;

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

std::vector<double> makeForkedDist(const Params& params) {
    return std::vector<double>(params.num_traits - 1, 0.5);
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
            } else {
                traitDist = makeConstrainedDist(params);
            }
            // loop over the traits, for each trait sample from the distribution. If the sample is smaller than the threshold, set the trait to 1
            for (size_t trait = 1; trait < params.num_traits; ++trait) {
                std::uniform_real_distribution<> dis(0, 1);
                if (dis(gen) < traitDist[trait - 1]) {
                    repertoires[agent][tree_idx][trait] = 1;
                } else break;
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
    // Adjusting for temperature in softmax
    double scaledSimilarity = expectedValues.similarity / params.temperature;
    double scaledProximal = expectedValues.proximal / params.temperature;
    
    double total = exp(scaledSimilarity) + exp(scaledProximal);
    double similarity_prob = exp(scaledSimilarity) / total;
    std::uniform_real_distribution<> dis(0, 1);
    
    if (dis(gen) < similarity_prob) {
        return Similarity;
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
    
    std::discrete_distribution<> distribution(traitWeights.begin(), traitWeights.end());

    size_t sampledTraitIndex = distribution(gen);
    return usefulTraits[sampledTraitIndex];
}


size_t Agents::learnSimilarity(size_t focalAgent, const std::vector<size_t>& usefulDemonstrators, size_t treeIndex, std::mt19937& gen) {
    std::vector<double> demonstratorWeights(usefulDemonstrators.size());
    std::vector<size_t> sampledTraits(usefulDemonstrators.size());
    double weightSum = 0.0;
    // if focal agent and demonstrator have the same trait, add to the weight
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
            demonstratorWeights[i] = std::inner_product(demoTraits.begin(), demoTraits.end(), repertoires[focalAgent][treeIndex].begin(), 0.0);
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

    // sample a trait from the weighted distribution
    std::discrete_distribution<> distribution(demonstratorWeights.begin(), demonstratorWeights.end());
    size_t sampledDemonstratorIndex = distribution(gen);
    return sampledTraits[sampledDemonstratorIndex];
}

size_t Agents::learnProximal(size_t focalAgent, const std::vector<size_t>& usefulDemonstrators, size_t treeIndex, std::mt19937& gen) {
    std::vector<double> demonstratorWeights(usefulDemonstrators.size());
    std::vector<size_t> sampledTraits(usefulDemonstrators.size());
    double weightSum = 0.0;

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
    
    // sample a trait from the weighted distribution
    std::discrete_distribution<> distribution(demonstratorWeights.begin(), demonstratorWeights.end());
    size_t sampledDemonstratorIndex = distribution(gen);
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
            feedback = payoffs[treeIndex][chosenTrait];
        } else {
            feedback = -1.0;
        }
        
        if (strategy == Similarity) {
            double predictionError = feedback - expectedValues[agentIndex].similarity;
            expectedValues[agentIndex].similarity += params.learning_rate * predictionError;
            chosenStrategy.push_back(Similarity);
        } else if (strategy == Proximal) {
            double predictionError = feedback - expectedValues[agentIndex].proximal;
            expectedValues[agentIndex].proximal += params.learning_rate * predictionError;
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
    if(debug >= 1) std::cout << "Update done" << '\n';
}

void Agents::learn(const Params& params, const std::vector<Tree>& trees) {
    std::random_device rd;
    std::mt19937 gen(rd());

    auto treeIndex = sampleIndex(params.num_trees, gen);
    auto agentIndex = sampleIndex(params.num_agents, gen);
    
    size_t chosenTrait;    
    std::uniform_real_distribution<> learn_dis(0, 1);
    Strategy strategy;

    auto unknownTraits = getUnknownTraits(agentIndex, treeIndex);
    if(debug >= 1) std::cout << "Unknown traits: " << unknownTraits.size() << '\n';
    if (unknownTraits.empty()) {
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
            if(debug >= 1) std::cout << "Strategy: " << strategy << '\n';
            auto demonstrators = sampleDemonstrators(agentIndex, params, gen);
            if(debug >= 1) std::cout << "Demonstrators: " << demonstrators.size() << '\n';
            if (strategy == Similarity) {
                // similarity learning
                if(debug >= 1) std::cout << "Similarity learning" << '\n';
                auto usefulDemonstrators = getUsefulDemonstrators(agentIndex, treeIndex, demonstrators);
                if(debug >= 1) std::cout << "Useful demonstrators: " << usefulDemonstrators.size() << '\n';
                
                if (usefulDemonstrators.empty()) {
                    strategy = None;
                    chosenTrait = 0;
                } else {
                    chosenTrait = learnSimilarity(agentIndex, usefulDemonstrators, treeIndex, gen);
                }
                if(debug >= 1) std::cout << "Similarity learning done" << '\n';
            } else {
                // proximal learning
                if(debug >= 1) std::cout << "Proximal learning" << '\n';
                auto usefulDemonstrators = getUsefulDemonstrators(agentIndex, treeIndex, demonstrators);
                if(debug >= 1) std::cout << "Useful demonstrators: " << usefulDemonstrators.size() << '\n';
                if (usefulDemonstrators.empty()) {
                    strategy = None;
                    chosenTrait = 0;
                } else {
                    chosenTrait = learnProximal(agentIndex, usefulDemonstrators, treeIndex, gen);
                }
                if(debug >= 1) std::cout << "Proximal learning done" << '\n';
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