#ifndef SNAPSHOT_HPP_INCLUDED_
#define SNAPSHOT_HPP_INCLUDED_

#include <sys/signal.h>
#include <vector>
#include <algorithm>
#include "mmapped_file.hpp"


#define ensure(cond)                                                          \
    do {                                                                      \
        if (!(cond)) {                                                        \
            fprintf(stderr, "Ensure `%s' in %s:%i (%s) failed.\n",            \
                    #cond, __FILE__, __LINE__, __FUNCTION__);                 \
            kill(0, SIGINT);                                                  \
        }                                                                     \
    } while (0)

#define ensure_valid_p(s, pidx)                                         \
    do                                                                  \
    {                                                                   \
        auto _pidx = pidx;                                              \
        ensure(_pidx >= 0 && static_cast<size_t>(_pidx) < (s).npart()); \
    } while (0)

//#if !defined(VERBOSE) || VERBOSE == 0
//#define verbose_printf(...) ((void) 0)
//#else
//#define verbose_printf printf
//#endif

typedef int particle_id;
typedef int bond_id;

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
        ensure(3 * npart() == pos.size());
        ensure(npart() + nproc() == boff.size());

        for (size_t i = 0; i < id.size(); ++i) {
            ensure_valid_p(*this, id[i]);
            ppermut[id[i]] = i;
        }
        ensure(std::count(std::begin(ppermut), std::end(ppermut), 0) == 1);
    }

    size_t npart() const { return id.size(); }
    size_t nproc() const { return pref.size(); }

    const double *pos_of_part(particle_id pid) const {
        ensure(pid >= 0 && static_cast<size_t>(pid) < npart());
        return &pos[3 * ppermut[pid]];
    }
    const double *vel_of_part(particle_id pid) const {
        ensure(pid >= 0 && static_cast<size_t>(pid) < npart());
        return &vel[3 * ppermut[pid]];
    }
};

struct Particle {
    particle_id id;
    const double *pos, *vel;
};

struct Bond {
    bond_id bid;
    particle_id pid;
    int npartners;
    const particle_id *partner_ids;
};

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
    ensure(static_cast<size_t>(glo_off) == s.bond.size());
}


inline double pdist(const snapshot& s, int pid1, int pid2)
{
    ensure(pid1 < s.npart() && pid2 < s.npart());

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

//double pdist1(const snapshot& s, int pidx1, int pid2)
//{
//    ensure_valid_p(s, pidx1); ensure_valid_p(s, pid2);
//
//    double pd = 0.0;
//    const auto pos1 = &s.pos[3 * pidx1];
//    auto pos2 = s.pos_of_part(pid2);
//
//    for (int d = 0; d < 3; ++d) {
//        auto dist = std::fabs(pos2[d] - pos1[d]);
//        if (dist > HALF_BOX_L)
//            dist = BOX_L - dist;
//        pd += dist * dist;
//    }
//
//    return std::sqrt(pd);
//}

inline double vec3len2(const double v[3])
{
    return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
}

inline double kinetic_energy(const snapshot& s, particle_id pid)
{
    ensure_valid_p(s, pid);
    return .5 * vec3len2(s.vel_of_part(pid));
}


#endif