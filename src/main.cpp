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
    for (size_t i = 0; i < params.num_iterations; ++i) {
        agents.learn(params, trees);
    }

    return Result{propConstrained, agents.computeProportion()};
}


int main() {
    Params params;
    
    std::vector<double> propConstrainedValues;
    for (double prop = 0.0; prop <= 1.0; prop += 0.05) {
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
    
    writeResultsToCsv(results, "../output.csv");
    std::cout << "Simulation complete" << '\n';
    return 0;
}