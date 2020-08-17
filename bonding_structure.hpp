// For license details see LICENSE.

#pragma once

#include <limits>
#include <map>
#include <vector>

#include "passert.hpp"
#include "types.hpp"
#include "uf.hpp"

/** Class storing relevant information about the bonding structure.
 * Slim wrapper around UnionFind.
 */
struct BondingStructure {

    BondingStructure(size_t npart) : uf(npart), npart(npart) {}

    /** Interface to register bond.
     */
    void add_bond(int i, int j) { uf.unionp(i, j); }

    /** Returns a vector containing all agglomerates. An agglomerate is, in
     * turn, a vector of particle ids.
     */
    std::vector<Agglomerate> agglomerates() {
        std::map<int, size_t> root_to_idx; //< Mapping of UF root node to index
                                           //in vector "agglomerates"
        std::vector<Agglomerate> agglomerates;

        p_assert(npart < std::numeric_limits<int>::max());
        for (int i = 0; i < static_cast<int>(npart); ++i) {
            auto r = uf.find(i);
            // We skip roots here in order to avoid size 1 clusters
            if (r == i)
                continue;
            auto [ii, unseen] = root_to_idx.emplace(r, agglomerates.size());
            if (unseen)
                agglomerates.push_back(Agglomerate{
                    r}); // Add the root here, since we skip it in the loop

            agglomerates[ii->second].push_back(i);
        }
        return agglomerates;
    }

    UnionFind uf;       //< UnionFind to compute the bonding structure
    const size_t npart; //< Number of particles
};
