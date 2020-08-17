// For license details see LICENSE.

#pragma once

#include <vector>
#include <map>
#include <algorithm>

/** Union find class with path pruning for integers.
 */
struct UnionFind {
    static constexpr const int root_marker = -1;
    std::vector<int> pred;

    UnionFind(size_t npart): pred(npart, root_marker) {}

    bool is_root(int i) {
        return pred[i] == root_marker;
    }

    int find(int i) {
        for (; !is_root(i); i = pred[i]);
        return i;
    }

    // Find and run function on every node on the path.
    // The function must have the following signature:
    // void f(int j) where j is the current node on the
    // path taken from i to its root node. The function
    // is allowed to modify pred[j], however, this does
    // not modify the path currently taken (the parent
    // node of j is determined before calling f).
    template <typename F>
    int findf(int i, F f) {
        for (; !is_root(i);) {
            // f(i) is allowed to change pred[i].
            auto j = pred[i];
            f(i);
            i = j;
        }
        return i;
    }


    void prune(int i, int root) {
        (void) findf(i, [root, this](int j){ pred[j] = root; });
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


/** Class storing relevant information about the bonding structure.
 */
struct BondingStructure {

    BondingStructure(size_t npart):
      uf(npart), npart(npart) {}

    /** Interface to register bond.
     */
    void add_bond(int i, int j) {
        uf.unionp(i, j);
    }

    typedef std::vector<int> agglo_type;

    /** Returns a vector containing all agglomerates. An agglonerate is, in turn,
     * a vector of particle ids.
     */
    std::vector<agglo_type> agglomerates() {
        std::map<int, size_t> root_to_idx; //< Mapping of UF root node to index in vector "agglomerates"
        std::vector<agglo_type> agglomerates;

        p_assert(npart < std::numeric_limits<int>::max());
        for (int i = 0; i < static_cast<int>(npart); ++i) {
            auto r = uf.find(i);
            // We skip roots here in order to avoid size 1 clusters
            if (r == i)
                continue;
            auto [ii, unseen] = root_to_idx.emplace(r, agglomerates.size());
            if (unseen)
                agglomerates.push_back(agglo_type{r}); // Add the root here, since we skip it in the loop

            agglomerates[ii->second].push_back(i);
        }
        return agglomerates;
    }

    UnionFind uf; //< UnionFind to compute the bonding structure
    const size_t npart; //< Number of particles
};

/** Parsed representation of a bond for storing it in memory.
 * Only used in case of "full bond storage".
 */
struct BondStore {
    BondStore(const BondReference &b) : bid(b.bid), pid(b.pid), npartners(b.npartners),
                               partner_ids(b.partner_ids, b.partner_ids + b.npartners)
    {}

    /** Returns the bond as non-owning "BondReference" data type.
     */
    operator BondReference() const {
        return {bid, pid, npartners, partner_ids.data()};
    }

    bond_id bid;
    particle_id pid;
    int npartners;
    std::vector<particle_id> partner_ids;
};

