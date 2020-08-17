// For license details see LICENSE.

#pragma once

#include <vector>

/** Union find class with path pruning for integers.
 */
struct UnionFind {
    static constexpr const int root_marker = -1;
    std::vector<int> pred;

    UnionFind(size_t npart) : pred(npart, root_marker) {}

    bool is_root(int i) { return pred[i] == root_marker; }

    int find(int i) {
        for (; !is_root(i); i = pred[i])
            ;
        return i;
    }

    // Find and run function on every node on the path.
    // The function must have the following signature:
    // void f(int j) where j is the current node on the
    // path taken from i to its root node. The function
    // is allowed to modify pred[j], however, this does
    // not modify the path currently taken (the parent
    // node of j is determined before calling f).
    template <typename F> int findf(int i, F f) {
        for (; !is_root(i);) {
            // f(i) is allowed to change pred[i].
            auto j = pred[i];
            f(i);
            i = j;
        }
        return i;
    }

    void prune(int i, int root) {
        (void)findf(i, [root, this](int j) { pred[j] = root; });
    }

    // Union and prune paths
    void unionp(int i, int j) {
        auto ri = find(i);
        auto rj = find(j);

        if (ri == rj)
            return;

        pred[rj] = ri;
        prune(j, ri);
    }
};
