#include "Utils.hpp"

#include <random>


size_t sampleIndex(const int& max, std::mt19937& gen) {
    std::uniform_int_distribution<> dis(0, max - 1);    
    return dis(gen);
}