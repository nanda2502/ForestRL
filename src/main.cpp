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
    }
    agents.printMeanEVs();

    result.timeSeriesData.reserve(agents.chosenStrategy.size());
    for (size_t i = 0; i < agents.chosenStrategy.size(); ++i) {
        TimeStepData timeStepData{
            propConstrained,
            i,  // Index in chosenStrategy is the timestep
            agents.chosenStrategy[i]
        };
        result.timeSeriesData.push_back(timeStepData);
    }

    return result;
}

int main(int argc, char* argv[]) {
    // optionally, accept input for showOutput
    int showOutput = 0;
    if (argc > 1) {
        showOutput = std::stoi(argv[1]);
    }

    Params params;
    
    std::vector<double> propConstrainedValues;
    for (double prop = 0.0; prop <= 0.0; prop += 0.25) {
        propConstrainedValues.push_back(prop);
    }
    
    std::vector<Result> results(propConstrainedValues.size());
    
    #pragma omp parallel for
    for (size_t i = 0; i < propConstrainedValues.size(); ++i) {
        double prop = propConstrainedValues[i];
        results[i] = runSimulation(params, prop, showOutput);
    }

    writeResultsToCsv(results, "../output.csv");
    std::cout << "Simulation complete" << '\n';
    return 0;
}