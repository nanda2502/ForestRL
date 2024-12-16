#include "Agents.hpp"
#include "Trees.hpp"
#include "Types.hpp"
#include "Utils.hpp"

#include <iostream>
#include <vector>

int debug = 0;

Result runSimulation(const Params& params, double propConstrained, int showOutput) {
    auto trees = initializeTrees(params, propConstrained);
    Agents agents(params, trees);
    Result result;
    
    for (size_t i = 0; i < params.num_iterations; ++i) {
        agents.learn(params, trees);
        if (i > 1e6 && showOutput == 1) debug = 1;
        if (i == 0) {
            agents.writeAveragesToCSV("repertoires_0.csv", params);
        }
        if (i == 2e6) {
            agents.writeAveragesToCSV("repertoires.csv", params);
        }
    }
    std::cout << "Proportion of constrained trees: " << propConstrained << '\n';
    agents.printMeanStratEVs();

    
    result.timeSeriesData.reserve(agents.chosenStrategy.size());
    for (size_t i = 0; i < agents.chosenStrategy.size(); ++i) {
        TimeStepData timeStepData{
            propConstrained,
            i,  // Index in chosenStrategy is the timestep
            agents.chosenStrategy[i],
            agents.ageAtTimestep[i],
            agents.chosenTree[i],
            agents.chosenAgent[i],
            agents.chosenTreeIndex[i],
            agents.receivedPayoffs[i]
        };
        result.timeSeriesData.push_back(timeStepData);
    }

    return result;
}

int main(int argc, char* argv[]) {
    int showOutput = 0;

    Params params = parseArgs(argc, argv);
    
    std::vector<double> propConstrainedValues;
    for (double prop = params.start_prop;
         prop <= params.end_prop;
         prop += params.prop_step) {
    propConstrainedValues.push_back(prop);
    }
    
    std::vector<Result> results(propConstrainedValues.size());
    
    #pragma omp parallel for
    for (size_t i = 0; i < propConstrainedValues.size(); ++i) {
        double prop = propConstrainedValues[i];
        results[i] = runSimulation(params, prop, showOutput);
    }

    std::string filename;
    if (params.update_trees) {
        filename = "./output_variable.csv";
    } else {
        filename = "./output_fixed.csv";
    }

    writeResultsToCsv(results, filename);
    std::cout << "Simulation complete" << '\n';
    return 0;
}