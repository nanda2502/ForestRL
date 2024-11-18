#include "Agents.hpp"
#include "Trees.hpp"
#include "Types.hpp"
#include "Utils.hpp"

#include <iostream>
#include <vector>
#include <execution>

int debug = 0;

Result runSimulation(const Params& params, double propConstrained) {
    auto trees = initializeTrees(params, propConstrained);
    Agents agents(params, trees);
    Result result;
    
    for (size_t i = 0; i < params.num_iterations; ++i) {
        agents.learn(params, trees);
        //if (i > 1e6) debug = 1;
    }

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

int main() {
    Params params;
    
    std::vector<double> propConstrainedValues;
    for (double prop = 0.0; prop <= 1.0; prop += 0.25) {
        propConstrainedValues.push_back(prop);
    }
    
    std::vector<Result> results(propConstrainedValues.size());
    
    std::transform(
        std::execution::par,
        propConstrainedValues.begin(),
        propConstrainedValues.end(),
        results.begin(),
        [&params](double prop) {
            return runSimulation(params, prop);
        }
    );

    #pragma omp parallel for
    for (size_t i = 0; i < propConstrainedValues.size(); ++i) {
        double prop = propConstrainedValues[i];
        results[i] = runSimulation(params, prop);
    }

    writeResultsToCsv(results, "../output.csv");
    std::cout << "Simulation complete" << '\n';
    return 0;
}