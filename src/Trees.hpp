#ifndef TREES_HPP
#define TREES_HPP

#include "Types.hpp"

#include <vector>
#include <random>

std::vector<Tree> initializeTrees(const Params& params, double propConstrained);

std::pair<std::vector<std::vector<double>>, double> generatePayoffs(const Params& params, std::mt19937& gen);

#endif // TREES_HPP