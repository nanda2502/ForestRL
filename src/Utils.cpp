#include "Utils.hpp"
#include "Types.hpp"

#include <random>
#include <iomanip>
#include <iostream>
#include <fstream>

size_t sampleIndex(const int& max, std::mt19937& gen) {
    std::uniform_int_distribution<> dis(0, max - 1);    
    return dis(gen);
}

void writeResultsToCsv(const std::vector<Result>& results, const std::string& filename) {
    std::ofstream outFile(filename);
    outFile << std::fixed << std::setprecision(3);
    outFile << "PropConstrained,Timestep,Strategy\n";
    
    for (const auto& result : results) {
        
        for (const auto& timeStep : result.timeSeriesData) {
            outFile << timeStep.propConstrained << "," 
                    << timeStep.timestep << ","
                    << timeStep.strategy << "\n";
        }
    }
}


