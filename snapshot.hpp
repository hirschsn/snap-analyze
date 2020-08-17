// For license details see LICENSE.

#pragma once

#include <algorithm>
#include <sys/signal.h>
#include <vector>

#include "mmapped_file.hpp"
#include "passert.hpp"
#include "span.hpp"
#include "types.hpp"

/** Ensure "pidx" is a valid particle given snapshot "s".
 */
#define ensure_valid_p(s, pidx)                                                \
    do {                                                                       \
        auto _pidx = pidx;                                                     \
        p_assert(_pidx >= 0 && static_cast<size_t>(_pidx) < (s).npart());      \
    } while (0)

template <typename Arr>
inline std::vector<int> create_inverse_permutation(const Arr &permut) {
    std::vector<int> inverse(permut.size(), 0);
    p_assert(permut.size() < std::numeric_limits<int>::max());
    for (int i = 0; i < static_cast<int>(permut.size()); ++i)
        inverse[permut[i]] = i;
    p_assert(std::count(std::begin(inverse), std::end(inverse), 0) == 1);
    return inverse;
}

/** Holds all data comprising a snapshot and defines accessor functions
 * to the data.
 * @see snapshot_iter()
 */
struct snapshot {
    const MFile<int> pref;
    const MFile<int> id;
    const MFile<double> pos, vel;
    const MFile<int> boff, bond;
    const MFile<int> bond_npartners;

    // Inverse permutation of "id"
    const std::vector<int> ppermut;

    const double box_l;
    const double half_box_l;

    snapshot(double box_l, std::string prefix)
        : pref((prefix + ".pref").c_str()), id((prefix + ".id").c_str()),
          pos((prefix + ".pos").c_str()), vel((prefix + ".vel").c_str()),
          boff((prefix + ".boff").c_str()), bond((prefix + ".bond").c_str()),
          bond_npartners(
              (prefix + ".head").c_str(),
              2 * sizeof(int)), // FIXME: sizeof(unsigned) + sizeof(size_t)
                                // Aktuell verwendetes Espresso schreibe noch
                                // unsigned + int raus
          ppermut(create_inverse_permutation(id)), box_l(box_l),
          half_box_l(box_l / 2.) {
        // Some sanity checks for the snapshot
        p_assert(3 * npart() == pos.size());
        p_assert(npart() + nproc() == boff.size());
        for (size_t i = 0; i < id.size(); ++i)
            ensure_valid_p(*this, id[i]);
    }

    size_t npart() const { return id.size(); }
    size_t nproc() const { return pref.size(); }

    const span3d pos_of_part(particle_id pid) const {
        p_assert(pid >= 0 && static_cast<size_t>(pid) < npart());
        return make_span3(&pos[3 * ppermut[pid]]);
    }
    const span3d vel_of_part(particle_id pid) const {
        p_assert(pid >= 0 && static_cast<size_t>(pid) < npart());
        return make_span3(&vel[3 * ppermut[pid]]);
    }
};

/** Function to iterate over all particles and bonds of a snapshot.
 * Use this template to parse a snapshot and provide the according
 * particle and bond callbacks.
 * particle callbacks must mave the signature: void(ParticleReference)
 * bond callbacks must have the signature: void(BondReference)
 * Where ParticleReference and BondReference are the struct types defined above.
 */
template <typename PCB, typename BCB>
void snapshot_iter(const snapshot &s, PCB particle_callback,
                   BCB bond_callback) {
    int glo_off = 0;
    for (size_t rank = 0; rank < s.nproc(); ++rank) {

        int pstart = s.pref[rank];
        int pend = rank == s.nproc() - 1 ? s.npart() : s.pref[rank + 1];
        int boff_start = s.pref[rank] + rank;
        int boff_end =
            rank + (rank == s.nproc() - 1 ? s.npart() : s.pref[rank + 1]);

        auto nlocalpart = pend - pstart;

        for (int i = 0; i < nlocalpart; ++i) {
            auto p = pstart + i;

            particle_callback(
                ParticleReference{s.id[p], &s.pos[3 * p], &s.vel[3 * p]});

            auto pb = boff_start + i;
            auto bond_start = glo_off + s.boff[pb];
            auto bond_end = glo_off + s.boff[pb + 1];
            for (int bidx = bond_start; bidx < bond_end;) {
                int bond_num = s.bond[bidx];
                bidx++;
                int npartners = s.bond_npartners[bond_num];
                bond_callback(
                    BondReference{bond_num, s.id[p], npartners, &s.bond[bidx]});
                bidx += npartners;
            }
        }
        glo_off += s.boff[boff_end];
    }

    // Check if snapshot is complete
    p_assert(static_cast<size_t>(glo_off) == s.bond.size());
}
