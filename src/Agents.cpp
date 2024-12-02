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
    expectedValues = std::vector<strategyExpectedValues>(params.num_agents, strategyExpectedValues{1.0, 1.0});
    lifetimes = std::vector<size_t>(params.num_agents);
    ages = std::vector<size_t>(params.num_agents, 0);
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
    return std::vector<double>(params.num_traits - 1, 0.8);
}


std::vector<double> makeFlatDist(const Params& params) {
    return std::vector<double>(params.num_traits - 1, 1.0/(params.num_traits - 1));
}

size_t sampleLifetime(const Params& params, std::mt19937& gen) {
    double mean =  params.lifetime_scale;

    std::poisson_distribution<size_t> lifetimeDist(mean);

    size_t baseLifetime = lifetimeDist(gen);

    return baseLifetime * params.num_trees;
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
                    };
                }
            } else {
                traitDist = makeConstrainedDist(params);
                for (size_t trait = 1; trait < params.num_traits; ++trait) {
                    std::uniform_real_distribution<> dis(0, 1);
                    if (dis(gen) < traitDist[trait - 1]) {
                        repertoires[agent][tree_idx][trait] = 1;
                    } else break; 
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

std::vector<size_t> Agents::sampleDemoTraits(size_t focalAgent, size_t treeIndex, const std::vector<size_t>& demonstrators, std::mt19937& gen) {
    std::vector<size_t> usefulTraits;
    const std::vector<size_t>& focalTraits = repertoires[focalAgent][treeIndex];

    for (size_t demonstrator : demonstrators) {
        const std::vector<size_t>& demoTraits = repertoires[demonstrator][treeIndex];

        std::vector<size_t> newTraits;
        for (size_t traitIndex = 0; traitIndex < demoTraits.size(); traitIndex++) {
            if (demoTraits[traitIndex] == 1 && focalTraits[traitIndex] == 0) {
                newTraits.emplace_back(traitIndex);
            }
        }
        if (!newTraits.empty()) {
            usefulTraits.emplace_back(newTraits[sampleIndex(newTraits.size(), gen)]);
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
    if (debug >=1) {
        std::vector<double> payoffsForTraits;
        payoffsForTraits.reserve(usefulTraits.size());
        for (size_t trait : usefulTraits) {
            payoffsForTraits.push_back(payoffs[treeIndex][trait]);
        }
        std::cout << "Payoffs for traits:\n";
        printVector(payoffsForTraits);
    }

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
    if (debug >= 1) std::cout << "Sampled trait: " << sampledTraits[sampledDemonstratorIndex] << '\n';
    return sampledTraits[sampledDemonstratorIndex];
}

int computeDelta(const std::vector<size_t>& focalTraits, const std::vector<size_t>& demoTraits) {
    return std::inner_product(focalTraits.begin(), focalTraits.end(), demoTraits.begin(), 0, 
    std::plus<>(),
    [](size_t focal_i, size_t demo_i) {
        return (focal_i == 1 && demo_i == 0) ? 1 : 0;
    }
    );
}

size_t Agents::learnProximalMarkov(size_t focalAgent, const std::vector<size_t>& usefulDemonstrators, size_t treeIndex, std::mt19937& gen) {
    auto focalTraits = repertoires[focalAgent][treeIndex];
    std::vector<double> traitWeights(focalTraits.size(), 0.0);
      
    
    for (size_t trait = 0; trait < focalTraits.size(); trait++) {
        if (focalTraits[trait] == 0) {
            for (size_t demonstrator : usefulDemonstrators) {
                const auto& demoTraits = repertoires[demonstrator][treeIndex];
                if (demoTraits[trait] == 1) {
                    auto delta = computeDelta(focalTraits, demoTraits);
                    if (delta > 0) {
                        traitWeights[trait] += std::pow(2.0, 1.0 - delta)/delta;
                    }
                }
            }
        }
    }

    if (debug >= 1) std::cout << "Trait weights:\n";
    if (debug >= 1) printVector(traitWeights);

    std::discrete_distribution<> distribution(traitWeights.begin(), traitWeights.end());
    size_t sampledTrait = distribution(gen);
    if (debug >= 1) std::cout << "Sampled trait: " << sampledTrait << '\n';
    return sampledTrait;
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
            feedback = payoffs[treeIndex][chosenTrait]; // payoffTotal;
        } else {
            feedback = 0.0;
        }
        if (debug >= 1) std::cout << "Feedback: " << feedback << '\n';
        if (strategy == Payoff) {
            double predictionError = feedback - expectedValues[agentIndex].payoff;
            expectedValues[agentIndex].payoff += params.learning_rate * predictionError;
            chosenStrategy.push_back(Payoff);
            ageAtTimestep.push_back(ages[agentIndex]);
        } else if (strategy == Proximal) {
            double predictionError = feedback - expectedValues[agentIndex].proximal;
            expectedValues[agentIndex].proximal += params.learning_rate * predictionError;
            chosenStrategy.push_back(Proximal);
            ageAtTimestep.push_back(ages[agentIndex]);
        } 
    }

    lifetimes[agentIndex] -= 1;
    ages[agentIndex] += 1;
    if (lifetimes[agentIndex] == 0) {
        if(debug >= 1) std::cout << "Agent reset" << '\n';
        //reset the agent's repertoire
        for (size_t i = 0; i < params.num_trees; ++i) {
            repertoires[agentIndex][i] = std::vector<size_t>(params.num_traits, 0);
            repertoires[agentIndex][i][0] = 1;
        }
        lifetimes[agentIndex] = sampleLifetime(params, gen);
        ages[agentIndex] = 0;
        expectedValues[agentIndex] = strategyExpectedValues{1.0, 1.0};
    }
}

void Agents::learn(const Params& params, const std::vector<Tree>& trees) {
    std::random_device rd;
    std::mt19937 gen(rd());

    auto agentIndex = sampleIndex(params.num_agents, gen);
    auto treeIndex = sampleIndex(params.num_trees, gen);

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
                auto usefulTraits = sampleDemoTraits(agentIndex, treeIndex, demonstrators, gen);

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
    if (debug >= 1) std::cout << "Agent " << agentIndex << " updated" << '\n';
}

double Agents::computeProportion() {
    int numElements = chosenStrategy.size() * 0.05;
    auto start = chosenStrategy.end() - numElements;
    double sum = std::accumulate(start, chosenStrategy.end(), 0);
    return sum / numElements;
}

void Agents::printMeanEVs() {
    size_t numAgents = expectedValues.size();

    double totalPayoffEV = 0.0;
    double totalProximalEV = 0.0;

    for (const auto& ev : expectedValues) {
        totalPayoffEV += ev.payoff;
        totalProximalEV += ev.proximal;
    }

    double meanPayoffEV = totalPayoffEV / numAgents;
    double meanProximalEV = totalProximalEV / numAgents;

    std::cout << "Mean Expected Value for Payoff Strategy: " << meanPayoffEV << '\n';
    std::cout << "Mean Expected Value for Proximal Strategy: " << meanProximalEV << '\n';
}


size_t Agents::sampleUnexploredTree(size_t agentIndex, std::mt19937& gen) {
    const auto& agentRepertoire = repertoires[agentIndex];
    std::vector<size_t> smallestTraitTrees;
    size_t minTraitCount = std::numeric_limits<size_t>::max();

    // Determine the minimum number of traits across all trees
    for (const auto& tree : agentRepertoire) {
        size_t traitCount = std::count(tree.begin(), tree.end(), 1);
        if (traitCount < minTraitCount) {
            minTraitCount = traitCount;
        }
    }

    // Collect all trees that have the minimum trait count
    for (size_t treeIndex = 0; treeIndex < agentRepertoire.size(); ++treeIndex) {
        const auto& tree = agentRepertoire[treeIndex];
        size_t traitCount = std::count(tree.begin(), tree.end(), 1);
        if (traitCount == minTraitCount) {
            smallestTraitTrees.push_back(treeIndex);
        }
    }

    // Sample randomly from the trees with the smallest number of traits
    std::uniform_int_distribution<size_t> treeDist(0, smallestTraitTrees.size() - 1);
    return smallestTraitTrees[treeDist(gen)];
}
