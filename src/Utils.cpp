#include "Utils.hpp"
#include "Types.hpp"

#include <random>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <unordered_map>

size_t sampleIndex(const int& max, std::mt19937& gen) {
    std::uniform_int_distribution<> dis(0, max - 1);    
    return dis(gen);
}

void writeResultsToCsv(const std::vector<Result>& results, const std::string& filename) {
    std::ofstream outFile(filename);
    outFile << std::fixed << std::setprecision(3);
    outFile << "PropConstrained,Timestep,Strategy,Age,Tree,Agent,TreeIndex,Payoff\n";
    
    for (const auto& result : results) {
        
        for (const auto& timeStep : result.timeSeriesData) {
            outFile << timeStep.propConstrained << "," 
                    << timeStep.timestep << ","
                    << timeStep.strategy << ","
                    << timeStep.age << ","
                    << timeStep.treeType << ","
                    << timeStep.agent << ","
                    << timeStep.treeIndex << ","
                    << timeStep.payoff << '\n';

        }
    }
}

Params parseArgs(int argc, char* argv[]) {
    Params params;
    std::unordered_map<std::string, std::string> argMap;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg.substr(0, 2) == "--" && i + 1 < argc) {
            std::string key = arg.substr(2);  // Remove leading "--"
            std::string value = argv[i + 1];
            argMap[key] = value;
            ++i;  // Skip the value for the next iteration
        }
    }

    try {
        if (argMap.find("num_agents") != argMap.end()) {
            params.num_agents = std::stoul(argMap["num_agents"]);
        }
        if (argMap.find("num_trees") != argMap.end()) {
            params.num_trees = std::stoul(argMap["num_trees"]);
        }
        if (argMap.find("num_demonstrators") != argMap.end()) {
            params.num_demonstrators = std::stoul(argMap["num_demonstrators"]);
        }
        if (argMap.find("num_traits") != argMap.end()) {
            params.num_traits = std::stoul(argMap["num_traits"]);
        }
        if (argMap.find("num_iterations") != argMap.end()) {
            params.num_iterations = std::stoul(argMap["num_iterations"]);
        }
        if (argMap.find("lifetime_scale") != argMap.end()) {
            params.lifetime_scale = std::stod(argMap["lifetime_scale"]);
        }
        if (argMap.find("learning_rate") != argMap.end()) {
            params.learning_rate = std::stod(argMap["learning_rate"]);
        }
        if (argMap.find("tree_learning_rate") != argMap.end()) {
            params.tree_learning_rate = std::stod(argMap["tree_learning_rate"]);
        }
        if (argMap.find("temperature") != argMap.end()) {
            params.temperature = std::stod(argMap["temperature"]);
        }
        if (argMap.find("tree_temperature") != argMap.end()) {
            params.tree_temperature = std::stod(argMap["tree_temperature"]);
        }
        if (argMap.find("innovation_rate") != argMap.end()) {
            params.innovation_rate = std::stod(argMap["innovation_rate"]);
        }
        if (argMap.find("constrained_payoff_scale") != argMap.end()) {
            params.constrained_payoff_scale = std::stod(argMap["constrained_payoff_scale"]);
        }

        if (argMap.find("start_prop") != argMap.end()) {
            params.start_prop = std::stod(argMap["start_prop"]);
        }

        if (argMap.find("end_prop") != argMap.end()) {
            params.end_prop = std::stod(argMap["end_prop"]);
        }

        if (argMap.find("prop_step") != argMap.end()) {
            params.prop_step = std::stod(argMap["prop_step"]);
        }

        if (argMap.find("update_trees") != argMap.end()) {
            params.update_trees = (argMap["update_trees"] == "true");
        }

    } catch (const std::exception& e) {
        std::cerr << "Error parsing arguments: " << e.what() << '\n';
    }

    return params;
}

