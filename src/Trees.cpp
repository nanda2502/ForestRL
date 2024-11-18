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

Tree makeForkedTree(const Params& params) {
    Tree tree(params.num_traits, std::vector<size_t>(params.num_traits, 0));

    if (params.num_traits > 1) {
        tree[0][1] = 1;
        if (params.num_traits > 2) {
            tree[0][2] = 1;
            for (size_t i = 3; i < params.num_traits; ++i) {
                if (i % 2 == 1) {
                    tree[i - 2][i] = 1;
                } else {
                    tree[i - 1][i] = 1;
                }
            }
        }
    }

    return tree;
}
Tree makeImbalancedTree(const Params& params) {
    Tree tree(params.num_traits, std::vector<size_t>(params.num_traits, 0));

    if (params.num_traits > 1) {
        // The number of direct children of the root is 2/3 of the total number of nodes excluding the root
        size_t total_nodes_excluding_root = params.num_traits - 1;
        size_t num_direct_children = (total_nodes_excluding_root * 2) / 3;
        
        // Connect the root to its direct children
        for (size_t i = 1; i <= num_direct_children; ++i) {
            tree[0][i] = 1;
        }

        // Form a single line with the remaining 1/3 of nodes
        if (params.num_traits > num_direct_children + 1) {
            size_t start_line_index = num_direct_children + 1;
            for (size_t i = start_line_index; i < params.num_traits; ++i) {
                tree[i - 1][i] = 1;
            }
        }
    }

    return tree;
}

std::vector<Tree> initializeTrees(const Params& params, double propConstrained) {
    auto numConstrainedTrees = static_cast<size_t>(params.num_trees * propConstrained);
    std::vector<Tree> trees(params.num_trees);

    for (size_t i = 0; i < numConstrainedTrees; ++i) {
        trees[i] = makeConstrainedTree(params);
    }

    // Forked is replacing flat here
    for (size_t i = numConstrainedTrees; i < params.num_trees; ++i) {
        trees[i] = makeForkedTree(params);
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

    // Generate params.num_trees random shuffles of these payoffs
    std::vector<std::vector<double>> payoffs(params.num_trees, std::vector<double>(params.num_traits));
    for (size_t tree = 0; tree < params.num_trees; ++tree) {
        std::vector<double> shuffled_payoffs = non_root_payoffs;
        std::shuffle(shuffled_payoffs.begin(), shuffled_payoffs.end(), gen);
        payoffs[tree][0] = 0.0;
        for (size_t i = 0; i < non_root_count; ++i) {
            payoffs[tree][i + 1] = shuffled_payoffs[i];
        }

    }
    return {payoffs, total};
}