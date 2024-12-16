#include "Trees.hpp"
#include <algorithm>

Tree makeConstrainedTree(const Params& params) {
    Tree tree(params.num_traits, std::vector<size_t>(params.num_traits, 0));

    for (size_t i = 0; i < params.num_traits -1; ++i) {
         tree[i][i+1] = 1;
    }

    return tree;
}

Tree makeFlatTree(const Params& params) {
    Tree tree(params.num_traits, std::vector<size_t>(params.num_traits, 0));

    for (size_t i = 1; i < params.num_traits; ++i) {
        tree[0][i] = 1;
    }

    return tree;
}

std::vector<Tree> initializeTrees(const Params& params, double propConstrained) {
    auto numConstrainedTrees = static_cast<size_t>(params.num_trees * propConstrained);
    std::vector<Tree> trees(params.num_trees);

    for (size_t i = 0; i < numConstrainedTrees; ++i) {
        trees[i] = makeConstrainedTree(params);
    }

    for (size_t i = numConstrainedTrees; i < params.num_trees; ++i) {
        trees[i] = makeFlatTree(params);
    }

    return trees;
}

std::pair<std::vector<std::vector<double>>, double> generatePayoffs(const Params& params, std::mt19937& gen) {
    constexpr double mean = 1.0;
   
    size_t non_root_count = params.num_traits - 1; 

    double spacing = (2.0 * mean) / (non_root_count + 1);
    std::vector<double> non_root_payoffs(non_root_count);
    double total = 0.0;
    // Generate equally spaced payoffs for non-root traits
    for (size_t i = 0; i < non_root_count; ++i) {
        non_root_payoffs[i] = spacing * (i + 1);
        total += non_root_payoffs[i];
    }

    std::vector<double> constrainedPayoffs(non_root_count);

    std::transform(non_root_payoffs.begin(), non_root_payoffs.end(), constrainedPayoffs.begin(), [&params](double payoff) {
        return payoff * params.constrained_payoff_scale;
    });

    // Payoffs for constrained trees
    std::vector<std::vector<double>> payoffs(params.num_trees, std::vector<double>(params.num_traits));

    for (size_t tree = 0; tree < (params.num_trees/2)+1; ++tree) {
        
        auto shuffledConstrainedPayoffs = constrainedPayoffs;
        std::shuffle(shuffledConstrainedPayoffs.begin(), shuffledConstrainedPayoffs.end(), gen);
        payoffs[tree][0] = 0.0;
        for (size_t i = 0; i < non_root_count; ++i) {
            payoffs[tree][i + 1] = shuffledConstrainedPayoffs[i];
        }
    }

    // Shuffle payoffs for flat trees
    for (size_t tree = (params.num_trees/2)+1; tree < params.num_trees; ++tree) {
        std::vector<double> shuffled_payoffs = non_root_payoffs;
        std::shuffle(shuffled_payoffs.begin(), shuffled_payoffs.end(), gen);
        payoffs[tree][0] = 0.0;
        for (size_t i = 0; i < non_root_count; ++i) {
            payoffs[tree][i + 1] = shuffled_payoffs[i];
        }
    }
    return {payoffs, total};
}