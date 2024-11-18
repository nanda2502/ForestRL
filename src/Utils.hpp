#ifndef UTILS_HPP
#define UTILS_HPP

#include "Types.hpp"

#include <random>
#include <iostream>

size_t sampleIndex(const int& max, std::mt19937& gen);
void writeResultsToCsv(const std::vector<Result>& results, const std::string& filename);

template <typename T>
void printVector(const std::vector<T>& v) {
    for (const auto& x : v) {
        std::cout << x << ' ';
    }
    std::cout << '\n';
}

#endif // UTILS_HPP