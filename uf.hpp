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
      uf(npart), npart(npart), max_bl_per_part(npart, 0.0),
      max_hbl_per_part(npart, 0.0), max_angle_per_part(npart, 0.0) {}

    /** Interface to register bond.
     */
    void add_bond(int i, int j, double dist) {
        uf.unionp(i, j);
        max_bl_per_part[i] = std::max(max_bl_per_part[i], dist);
        max_bl_per_part[j] = std::max(max_bl_per_part[j], dist);
    }

    /** Interface to register bond length of a harmonic bond.
     * If you call this, also call add_bond().
     */
    void add_hbond(int i, int j, double dist) {
        max_hbl_per_part[i] = std::max(max_hbl_per_part[i], dist);
        max_hbl_per_part[j] = std::max(max_hbl_per_part[j], dist);
    }

    /** Interface to register angles of an angular bond.
     * If you call this, also call add_bond().
     */
    void add_angle(int i, int j, int k, double angle) {
        max_angle_per_part[i] = std::max(max_angle_per_part[i], std::fabs(angle));
        max_angle_per_part[j] = std::max(max_angle_per_part[j], std::fabs(angle));
        max_angle_per_part[k] = std::max(max_angle_per_part[k], std::fabs(angle));
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

    /** Maximal bond lengths, maximal harmonic bond lengths and maximal
     * angle of angle bonds per particle.
     */
    std::vector<double> max_bl_per_part, max_hbl_per_part, max_angle_per_part;
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

typedef std::vector<std::vector<BondStore>> FullBondStorage;
