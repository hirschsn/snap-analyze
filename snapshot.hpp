// For license details see LICENSE.

#pragma once

#include <sys/signal.h>
#include <vector>
#include <algorithm>
#include "mmapped_file.hpp"


[[noreturn]] void __passert_fail(const char *expr, const char *file, int line, const char *function)
{
    std::fprintf(stderr, "p_assert assertion failed: `%s' in %s:%i (%s)", expr, file, line, function);
    std::abort();
}

/** Assert-equivalent that is *not* a no-op if NDEBUG is set.
 */
#define p_assert(cond)                                                          \
    ((cond)? (void) 0: __passert_fail(#cond, __FILE__, __LINE__, __FUNCTION__))

/** Ensure "pidx" is a valid particle given snapshot "s".
 */
#define ensure_valid_p(s, pidx)                                         \
    do                                                                  \
    {                                                                   \
        auto _pidx = pidx;                                              \
        p_assert(_pidx >= 0 && static_cast<size_t>(_pidx) < (s).npart()); \
    } while (0)

typedef int particle_id;
typedef int bond_id;

/** Holds all data comprising a snapshot and defines accessor functions
 * to the data.
 * @see snapshot_iter()
 */
struct snapshot {
    MFile<int> pref;
    MFile<int> id;
    MFile<double> pos, vel;
    MFile<int> boff, bond;
    MFile<int> bond_npartners;

    // Inverse permutation of "id"
    std::vector<int> ppermut;


    snapshot(std::string prefix): pref((prefix + ".pref").c_str()),
                                  id((prefix + ".id").c_str()),
                                  pos((prefix + ".pos").c_str()),
                                  vel((prefix + ".vel").c_str()),
                                  boff((prefix + ".boff").c_str()),
                                  bond((prefix + ".bond").c_str()),
                                  bond_npartners((prefix + ".head").c_str(),
                                                 2 * sizeof(int)), // FIXME: sizeof(unsigned) + sizeof(size_t)
                                                                   // Aktuell verwendetes Espresso schreibe noch unsigned + int raus
                                  ppermut(id.size(), 0) {
        p_assert(3 * npart() == pos.size());
        p_assert(npart() + nproc() == boff.size());

        for (size_t i = 0; i < id.size(); ++i) {
            ensure_valid_p(*this, id[i]);
            ppermut[id[i]] = i;
        }
        p_assert(std::count(std::begin(ppermut), std::end(ppermut), 0) == 1);
    }

    size_t npart() const { return id.size(); }
    size_t nproc() const { return pref.size(); }

    const double *pos_of_part(particle_id pid) const {
        p_assert(pid >= 0 && static_cast<size_t>(pid) < npart());
        return &pos[3 * ppermut[pid]];
    }
    const double *vel_of_part(particle_id pid) const {
        p_assert(pid >= 0 && static_cast<size_t>(pid) < npart());
        return &vel[3 * ppermut[pid]];
    }
};

/** Non-owning interface for access to particles.
 */
struct Particle {
    particle_id id;
    const double *pos, *vel;
};

/** Non-owning interface for access to bonds.
 */
struct Bond {
    bond_id bid;
    particle_id pid;
    int npartners;
    const particle_id *partner_ids;
};

/** Function to iterate over all particles and bonds of a snapshot.
 * Use this template to parse a snapshot and provide the according
 * particle and bond callbacks.
 * particle callbacks must mave the signature: void(Particle)
 * bond callbacks must have the signature: void(Bond)
 * Where Particle and Bond are the struct types defined above.
 */
template <typename PCB, typename BCB>
void snapshot_iter(const snapshot& s, PCB particle_callback, BCB bond_callback)
{
    int glo_off = 0;
    for (size_t rank = 0; rank < s.nproc(); ++rank) {

        int pstart = s.pref[rank];
        int pend = rank == s.nproc() - 1? s.npart(): s.pref[rank + 1];
        int boff_start = s.pref[rank] + rank;
        int boff_end = rank + (rank == s.nproc() - 1
                        ? s.npart()
                        : s.pref[rank + 1]);

        auto nlocalpart = pend - pstart;

        for (int i = 0; i < nlocalpart; ++i) {
            auto p = pstart + i;

            particle_callback(Particle{s.id[p], &s.pos[3 * p], &s.vel[3 * p]});

            auto pb = boff_start + i;
            auto bond_start = glo_off + s.boff[pb];
            auto bond_end = glo_off + s.boff[pb + 1];
            for (int bidx = bond_start; bidx < bond_end; ) {
                int bond_num = s.bond[bidx];
                bidx++;
                int npartners = s.bond_npartners[bond_num];
                bond_callback(Bond{bond_num, s.id[p], npartners, &s.bond[bidx]});
                bidx += npartners;
            }
        }
        glo_off += s.boff[boff_end];
    }

    // Check if snapshot is complete
    p_assert(static_cast<size_t>(glo_off) == s.bond.size());
}


/** Calculate the distance between two particles given their ids.
 */
inline double pdist(const snapshot& s, int pid1, int pid2)
{
    p_assert(pid1 < s.npart() && pid2 < s.npart());

    double pd = 0.0;
    auto pos1 = s.pos_of_part(pid1);
    auto pos2 = s.pos_of_part(pid2);

    for (int d = 0; d < 3; ++d) {
        auto dist = std::fabs(pos2[d] - pos1[d]);
        if (dist > HALF_BOX_L)
            dist = BOX_L - dist;
        pd += dist * dist;

    }

    return std::sqrt(pd);
}

/** Squared vector norm in 3d
 */
inline double vec3len2(const double v[3])
{
    return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
}

/** Calculates the kinetic energy of a particle given its id.
 */
inline double kinetic_energy(const snapshot& s, particle_id pid)
{
    ensure_valid_p(s, pid);
    return .5 * vec3len2(s.vel_of_part(pid));
}
