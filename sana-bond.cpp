// For license details see LICENSE.
 
#include <iostream>
#include <cmath>
#include <utility>
#include <map>
#include <iterator>
#include <cmath>
#include <numeric>
#include <string>
#include <optional>
#include "box.hpp"
#include "stat.hpp"
#include "snapshot.hpp"
#include "df.hpp"
#include "uf.hpp"

bool is_broken_angle(double angle)
{
    return angle <= -.3 || angle >= .3;
}

std::vector<Vec3d> ids_to_poss(const snapshot &s, const std::vector<int> &agg)
{
    std::vector<Vec3d> poss(agg.size());

    for (size_t i = 0; i < agg.size(); ++i) {
        for (int d = 0; d < 3; ++d)
            poss[i][d] = s.pos[3 * s.ppermut[agg[i]] + d];
    }

    return poss;
}
std::vector<Vec3d> ids_to_vels(const snapshot &s, const std::vector<int> &agg)
{
    std::vector<Vec3d> vels(agg.size());

    for (size_t i = 0; i < agg.size(); ++i) {
        for (int d = 0; d < 3; ++d)
            vels[i][d] = s.vel[3 * s.ppermut[agg[i]] + d];
    }

    return vels;
}

std::map<int, int> size_histogram(const std::vector<std::vector<int>>& aggs)
{
    std::map<int, int> hist;
    for (const auto& agg: aggs) {
        auto [ii, unseen] = hist.emplace(agg.size(), 0);
        ii->second++;
    }

    return hist;
}

void print_agglomerate_sizes_histogram(const std::vector<std::vector<int>>& aggs)
{
    auto hist = size_histogram(aggs);

    for (auto [size, n]: hist) {
        std::cout << "Size " << size << ": " << n;
        std::cout << std::endl;
    }
}

void print_agglomerate_sizes(const std::vector<std::vector<int>>& aggs)
{
    for (const auto& agg: aggs) {
        std::cout << agg.size() << std::endl;
    }
}

void print_agglomerate_sizes_with_maxdists(const snapshot& s, const std::vector<std::vector<int>>& aggs)
{
    for (const auto& agg: aggs) {
        auto pos = ids_to_poss(s, agg);
        auto m = maxdist(pos);
        std::cout << agg.size() << " " << m << std::endl;
    }
}

void print_agglomerate_sizes_with_maxbl(const BondingStructure& bs, const std::vector<std::vector<int>>& aggs)
{
    for (const auto& agg: aggs) {
        double ml = 0.0;
        for (int id: agg)
            ml = std::max(ml, bs.max_bl_per_part[id]);
        std::cout << agg.size() << " " << ml << std::endl;
    }
}

void print_agglomerate_sizes_with_maxangle(const BondingStructure& bs, const std::vector<std::vector<int>>& aggs)
{
    for (const auto& agg: aggs) {
        double ma = 0.0;
        for (int id: agg)
            ma = std::max(ma, bs.max_angle_per_part[id]);
        std::cout << agg.size() << " " << ma << std::endl;
    }
}

void print_agglomerate_sizes_with_all(const snapshot& s, const BondingStructure& bs, const std::vector<std::vector<int>>& aggs)
{
    for (const auto& agg: aggs) {
        auto pos = ids_to_poss(s, agg);
        double ma = 0.0, ml = 0.0;
        for (int id: agg) {
            ma = std::max(ma, bs.max_angle_per_part[id]);
            ml = std::max(ml, bs.max_bl_per_part[id]);
        }
        auto md = maxdist(pos);
        std::cout << agg.size() << " " << md << " " << ml << " " << ma << std::endl;
    }
}

void print_agglomerate_of_size(const snapshot& s, const std::vector<std::vector<int>>& aggs, size_t size, bool print_pos, bool print_vel, bool print_id=false)
{
    int i = 0;
    for (const auto& agg: aggs) {
        if (agg.size() == size) {
            printf("Agglomerate %i of size %zu\n", i++, size);
            std::vector<Vec3d> pos, vel;

            if (print_pos)
                pos = ids_to_poss(s, agg);
            if (print_vel)
                vel = ids_to_vels(s, agg);
            for (int i = 0; i < agg.size(); ++i) {
                if (print_id)
                    printf("%i ", agg[i]);
                if (print_pos)
                    printf("%lf %lf %lf", pos[i][0], pos[i][1], pos[i][2]);
                if (print_pos && print_vel)
                    printf(" ");
                if (print_vel)
                    printf("%lf %lf %lf", vel[i][0], vel[i][1], vel[i][2]);
                if (print_pos || print_vel)
                    printf("\n");
            }
        }
    }
}

void print_agglomerate_tcl(const snapshot& s, const std::vector<std::vector<int>>& aggs, size_t size, bool print_pos, bool print_vel)
{
    int i = 0;
    for (const auto& agg: aggs) {
        if (agg.size() == size) {
            std::vector<Vec3d> pos, vel;

            if (print_pos)
                pos = ids_to_poss(s, agg);
            if (print_vel)
                vel = ids_to_vels(s, agg);
            for (int i = 0; i < agg.size(); ++i) {
                printf("\t{%i", i);
                if (print_pos)
                    printf(" %lf %lf %lf", pos[i][0], pos[i][1], pos[i][2]);
                if (print_vel)
                    printf(" %lf %lf %lf", vel[i][0], vel[i][1], vel[i][2]);
                printf("}\n");
            }
        }
    }
}

void print_bonding_info_tcl(const FullBondStorage& fbs, const std::vector<int>& agg, const std::map<int, int>& id_rename)
{
    auto r = [&id_rename](int i){return id_rename.at(i); };

    for (auto pid: agg) {
        const auto& bl = fbs[pid];
        printf("\t{%i { ", r(pid));
        for (const auto& b: bl) {
            if (b.npartners == 1)
                printf("{%i %i} ", b.bid, r(b.partner_ids[0]));
            else if (b.npartners == 2)
                printf("{%i %i %i} ", b.bid, r(b.partner_ids[0]), r(b.partner_ids[1]));
            else
                throw std::runtime_error("b.nparters: " + std::to_string(b.npartners));
        }
        printf("} }\n");
    }
}

void print_agglomerate_raw(const snapshot& s, const FullBondStorage& fbs, const std::vector<std::vector<int>>& aggs, size_t size)
{
    // ID renaming - make them sequential starting from 0
    std::map<int, int> id_rename;
    for (const auto& agg: aggs) {
        if (agg.size() == size) {
            int idno = 0;
            for (auto id: agg) {
                id_rename[id] = idno++;
            }
        }
    }
    printf("{particles {id pos v}\n");
    print_agglomerate_tcl(s, aggs, size, true, true);
    printf("}\n");
    printf("{bonds\n");
    for (const auto& agg: aggs) {
        if (agg.size() == size) {
            print_bonding_info_tcl(fbs, agg, id_rename);
        }
    }
    printf("}\n");
}

void print_agglomerates_to_files(const snapshot& s, const BondingStructure& bs, const std::vector<std::vector<int>>& aggs)
{
    std::map<int, int> num;
    for (const auto& agg: aggs) {
        auto [ii, unseen] = num.emplace(agg.size(), 0);
        ii->second++;
        auto pos = ids_to_poss(s, agg);
        auto fn = std::to_string(agg.size()) + "_POS_" + std::to_string(ii->second);
        FILE *f = fopen(fn.c_str(), "w");

        for (int i = 0; i < agg.size(); ++i) {
            fprintf(f, "%lf %lf %lf\n", pos[i][0], pos[i][1], pos[i][2]);
        }
        fclose(f);
    }
}

void print_agglomerates_of_sizes(const snapshot& s, const BondingStructure& bs, const std::vector<std::vector<int>>& aggs, const std::vector<size_t>& wanted)
{
    std::map<int, int> num;

    for (const auto& agg: aggs) {
        if (std::find(std::begin(wanted), std::end(wanted), agg.size()) != std::end(wanted)) {
            auto [ii, unseen] = num.emplace(agg.size(), 0);
            ii->second++;
            auto pos = ids_to_poss(s, agg);
            auto vel = ids_to_vels(s, agg);


            auto fn = std::to_string(agg.size()) + "_POSVEL_" + std::to_string(ii->second);
            auto fn_bhb = fn + "_BHB";
            auto fn_bab = fn + "_BAB";
            auto fn_bbb = fn + "_BBB";
            FILE *f = fopen(fn.c_str(), "w");
            FILE *fh = fopen(fn_bhb.c_str(), "w");
            FILE *fa = fopen(fn_bab.c_str(), "w");
            FILE *fb = fopen(fn_bbb.c_str(), "w");

            for (int i = 0; i < agg.size(); ++i) {
                fprintf(f, "%lf %lf %lf %lf %lf %lf\n", pos[i][0], pos[i][1], pos[i][2], vel[i][0], vel[i][1], vel[i][2]);

                auto bh = bs.max_hbl_per_part[agg[i]] >= 2.0;
                auto ba = is_broken_angle(bs.max_angle_per_part[agg[i]]);
                if (bh)
                    fprintf(fh, "%lf %lf %lf %lf %lf %lf\n", pos[i][0], pos[i][1], pos[i][2], vel[i][0], vel[i][1], vel[i][2]);
                if (ba)
                    fprintf(fa, "%lf %lf %lf %lf %lf %lf\n", pos[i][0], pos[i][1], pos[i][2], vel[i][0], vel[i][1], vel[i][2]);
                if (bh && ba)
                    fprintf(fb, "%lf %lf %lf %lf %lf %lf\n", pos[i][0], pos[i][1], pos[i][2], vel[i][0], vel[i][1], vel[i][2]);
            }
            fclose(f);
            fclose(fh);
            fclose(fa);
        }
    }
}

void print_agglomerate_of_size_broken_bonds(const BondingStructure& bs, const snapshot& s, const std::vector<std::vector<int>>& aggs, size_t size)
{
    int i = 0;
    for (const auto& agg: aggs) {
        if (agg.size() == size) {
            auto pos = ids_to_poss(s, agg);
            printf("Agglomerate %i of size %zu\n", i++, size);
            for (size_t i = 0; i < agg.size(); ++i) {
                if (bs.max_bl_per_part[agg[i]] > 2.0) {
                    printf("%lf %lf %lf\n", pos[i][0], pos[i][1], pos[i][2]);
                }
            }
        }
    }
}

void print_agglomerate_of_size_broken_angle_bonds(const BondingStructure& bs, const snapshot& s, const std::vector<std::vector<int>>& aggs, size_t size)
{
    int i = 0;
    for (const auto& agg: aggs) {
        if (agg.size() == size) {
            auto pos = ids_to_poss(s, agg);
            printf("Agglomerate %i of size %zu\n", i++, size);
            for (size_t i = 0; i < agg.size(); ++i) {
                if (is_broken_angle(bs.max_angle_per_part[agg[i]])) {
                    printf("%lf %lf %lf\n", pos[i][0], pos[i][1], pos[i][2]);
                }
            }
        }
    }
}

void print_broken_angle_bonds_with_positions(const BondingStructure& bs, const snapshot& s, const std::vector<std::vector<int>>& aggs)
{
    for (const auto& agg: aggs) {

        double ma = 0.0;
        int pid = -1;
        for (int id: agg) {
            if (bs.max_angle_per_part[id] > ma) {
                ma = bs.max_angle_per_part[id];
                pid = id;
            }
        }
        if (ma < 0.3)
            continue;
        printf("Agglomerate of size %zu with max angle: %lf at particle %i\n", agg.size(), ma, pid);
        auto pos = ids_to_poss(s, agg);
        for (size_t i = 0; i < agg.size(); ++i)
            printf("%i %lf %lf %lf\n", agg[i], pos[i][0], pos[i][1], pos[i][2]);
    }
}

void print_all_broken_bonds(const BondingStructure& bs, const snapshot& s, const std::vector<std::vector<int>>& aggs)
{
    for (const auto& agg: aggs) {
        double ma = 0.0, ml = 0.0;
        for (int id: agg) {
            ma = std::max(ma, bs.max_angle_per_part[id]);
            ml = std::max(ml, bs.max_hbl_per_part[id]);
        }

        if (is_broken_angle(ma) && ml < 2.0)
            continue;


        std::vector<int> blp, bap, vint;
        for (int id: agg) {
            bool bl = bs.max_hbl_per_part[id] >= 2.0;
            bool ba = is_broken_angle(bs.max_angle_per_part[id]);
            if (bl)
                blp.push_back(id);
            if (ba)
                bap.push_back(id);
        }

        std::sort(blp.begin(), blp.end());
        std::sort(bap.begin(), bap.end());
        std::set_intersection(blp.begin(), blp.end(), bap.begin(), bap.end(), std::back_inserter(vint));

        std::cout << agg.size() << " " << ml << " " << ma << " broken-by-bl: " << blp.size() << " broken-by-phi: " << bap.size() << " in-both-sets: " << vint.size() << std::endl;
    }
}


void print_all_broken_bonds_verbose(const BondingStructure& bs, const snapshot& s, const std::vector<std::vector<int>>& aggs)
{
    for (const auto& agg: aggs) {
        double ma = 0.0, ml = 0.0;
        for (int id: agg) {
            ma = std::max(ma, bs.max_angle_per_part[id]);
            ml = std::max(ml, bs.max_hbl_per_part[id]);
        }

        if (!is_broken_angle(ma) && ml < 2.0)
            continue;

        std::cout << agg.size() << " " << ml << " " << ma << " broken-by-bl:";
        for (int id: agg) {
            if (bs.max_hbl_per_part[id] >= 2.0)
                std::cout << " " << id << " (" << bs.max_hbl_per_part[id] << ")";
        }

        std::cout << " broken-by-angle:";
        for (int id: agg) {
            if (is_broken_angle(bs.max_angle_per_part[id]))
                std::cout << " " << id << " (" << bs.max_angle_per_part[id] << ")";
        }
        std::cout << std::endl;
    }
}

void print_agglomerate_of_pid(const snapshot& s, const std::vector<std::vector<int>>& aggs, int pid)
{
  for (const auto& agg: aggs) {
    if (std::find(std::begin(agg), std::end(agg), pid) == std::end(agg))
      continue;

    auto poss = ids_to_poss(s, agg);
    for (size_t i = 0; i < agg.size(); ++i) {
      std::cout << agg[i] << ": " << poss[i][0] << " " << poss[i][1] << " " << poss[i][2] << "\n";
    }
    
    return;
  }
}



void print_agglomerate_dfs(const snapshot& s, const std::vector<std::vector<int>>& aggs)
{
//    ScalarStatistics<double> dfstat, maxdists;
//
    auto calc_df_radog = [&s](const std::vector<int> &agg){
        return calc_df(ids_to_poss(s, agg));
    };

//#pragma omp parallel for schedule(guided)
//    for (size_t i = 0; i < aggs.size(); ++i) {
//        const auto &agg = aggs[i];
//
//        if (agg.size() >= 15) {
//            auto [md, df_radog] = calc_df_radog(agg);
//#pragma omp critical
//            {
//                dfstat.sample(df_radog);
//                maxdists.sample(md);
//            }
//        }
//    }
    for (const auto &agg: aggs) {
        if (agg.size() < 15)
            continue;
        auto [radog, df_radog] = calc_df_radog(agg);
        //dfstat.sample(df_radog);

        std::cout << agg.size() << " " << radog << " " << df_radog << std::endl;
    }
}


double dround(double d)
{
    return std::floor(d + .5);
}

template <typename T>
Vec3d get_mi_vector(const T& a, const T& b) {
  Vec3d res;
  for (int i = 0; i < 3; i++) {
    res[i] = a[i] - b[i];
    if (std::fabs(res[i]) > HALF_BOX_L)
      res[i] -= dround(res[i] / BOX_L) * BOX_L;
  }

  return res;
}

/** calculates the squared length of a vector */
template <typename T>
double vlen(T const &v) {
  double d2 = 0.0;
  int i;
  for (i = 0; i < 3; i++)
    d2 += v[i] * v[i];
  return std::sqrt(d2);
}

double calc_bond_angle(const snapshot& s, const Bond& b)
{
    auto pos_mid = s.pos_of_part(b.pid);
    auto pos_left = s.pos_of_part(b.partner_ids[0]);
    auto pos_right = s.pos_of_part(b.partner_ids[1]);

    /* vector from p_left to p_mid */
    auto vec1 = get_mi_vector(pos_mid, pos_left);
    auto d1i = 1.0 / vlen(vec1);
    for (int j = 0; j < 3; j++)
        vec1[j] *= d1i;

    /* vector from p_mid to p_right */
    auto vec2 = get_mi_vector(pos_right, pos_mid);
    auto d2i = 1.0 / vlen(vec2);
    for (int j = 0; j < 3; j++)
        vec2[j] *= d2i;

    /* scalar product of vec1 and vec2 */
    auto cosine = std::inner_product(std::begin(vec1), std::end(vec1), std::begin(vec2), 0.0);
    
    constexpr const double bend = 1000.0;
    constexpr const double PI = 3.141592653589793;
    const double phi0 = b.bid * PI / 180.0;
    constexpr const double TINY_COS_VALUE = 0.9999999999;
    constexpr const double TINY_SIN_VALUE = 1e-10;

    double fac = bend;

    if (cosine > TINY_COS_VALUE)
        cosine = TINY_COS_VALUE;
    if (cosine < -TINY_COS_VALUE)
        cosine = -TINY_COS_VALUE;
    double phi = acos(-cosine);


    //std::cout << "Bond with bond_id " << b.bid << " and phi0 " << phi0 << " has angle phi = " << phi << std::endl;

    return phi - phi0;

    //double sinphi = sin(phi);
    //if (sinphi < TINY_SIN_VALUE)
    //    sinphi = TINY_SIN_VALUE;
    //fac *= (phi - phi0) / sinphi;

    //for (j = 0; j < 3; j++)
    //{
    //    double f1 = fac * (cosine * vec1[j] - vec2[j]) * d1i;
    //    double f2 = fac * (cosine * vec2[j] - vec1[j]) * d2i;

    //    force1[j] = (f1 - f2);
    //    force2[j] = -f1;
    //}

    //return 0;
}

std::optional<BondStore> find_pair_bond(const FullBondStorage& fbs, particle_id pid1, particle_id pid2)
{
    const auto& bl1 = fbs[pid1];
    const auto& bl2 = fbs[pid2];

    auto it = std::find_if(std::begin(bl1), std::end(bl1),
                           [pid2](const BondStore& bs){
                               return bs.npartners == 1 && bs.partner_ids[0] == pid2;
                           });
    
    if (it != std::end(bl1)) {
        return *it;
    }

    it = std::find_if(std::begin(bl2), std::end(bl2),
                           [pid1](const BondStore& bs){
                               return bs.npartners == 1 && bs.partner_ids[0] == pid1;
                           });
    
    if (it != std::end(bl2)) {
        return *it;
    }

    return {};
}

std::optional<BondStore> find_angle_bond(const FullBondStorage& fbs, particle_id pid1, particle_id pid2, particle_id pid3)
{
    const auto& bl1 = fbs[pid1];

    auto it = std::find_if(std::begin(bl1), std::end(bl1),
                           [pid2, pid3](const BondStore& bs){
                               return bs.npartners == 2
                                      && ((bs.partner_ids[0] == pid2 && bs.partner_ids[1] == pid3)
                                          || (bs.partner_ids[0] == pid3 && bs.partner_ids[1] == pid2));
                           });
    
    if (it != std::end(bl1)) {
        return *it;
    }

    return {};
}

void check_angle_bonds(const BondingStructure &bs, const snapshot &s, const std::vector<std::vector<int>> &aggs, const FullBondStorage &fbs) {
    for (const auto& agg: aggs) {
        for (int pid: agg) {
            for (const auto& b: fbs[pid]) {
                // Angular bond
                if (b.npartners == 2) {
                    auto bb = static_cast<Bond>(b);
                    auto angle = calc_bond_angle(s, b);

                    if (is_broken_angle(angle)) {
                        auto pid1 = bb.partner_ids[0];
                        auto pid2 = bb.partner_ids[1];
                        printf("Broken angle bond (%i) at particle %i angle %lf partners: %i %i", bb.bid, bb.pid, angle, pid1, pid2);

                        printf(" pair bonds:");
                        // Search for the opposite and neighboring pair bonds
                        std::array<std::pair<std::optional<BondStore>, double>, 3> obs = {{
                            { find_pair_bond(fbs, pid1, pid2), pdist(s, pid1, pid2) },
                            { find_pair_bond(fbs, pid, pid1), pdist(s, pid, pid1) },
                            { find_pair_bond(fbs, pid, pid2), pdist(s, pid, pid2) },
                        }};

                        for (const auto [ob, dist]: obs) {
                            if (ob) {
                                if (dist >= 2.0)
                                    printf(" broken (%i: %lf)", ob->bid, dist);
                                else
                                    printf(" valid (%i: %lf)", ob->bid, dist);
                            } else {
                                    printf(" none ()");
                            }
                        }

                        printf(" angle bonds:");
                        std::array<std::array<particle_id, 3>, 2> triangles = {{
                            {{ pid1, pid, pid2 }},
                            {{ pid2, pid, pid1 }},
                        }};

                        for (const auto [t1, t2, t3]: triangles) {
                            auto ob = find_angle_bond(fbs, t1, t2, t3);
                            if (ob) {
                                auto angle = calc_bond_angle(s, *ob);
                                if (is_broken_angle(angle))
                                    printf(" broken (%i: %lf)", ob->bid, angle);
                                else
                                    printf(" valid (%i: %lf)", ob->bid, angle);
                            } else {
                                    printf(" none ()");
                            }
                        }

                        putchar('\n');
                    }
                }
            }
        }
    }
}
void track_particle(const BondingStructure &bs, const snapshot &s, const std::vector<std::vector<int>> &aggs, const FullBondStorage &fbs, int wanted_id) {
    for (const auto& agg: aggs) {
        for (int pid: agg) {
            if (pid != wanted_id)
                continue;

            for (const auto& b: fbs[pid]) {
                // Angular bond
                if (b.npartners == 2) {
                    auto bb = static_cast<Bond>(b);
                    auto angle = calc_bond_angle(s, b);

                    auto pid1 = bb.partner_ids[0];
                    auto pid2 = bb.partner_ids[1];
                    printf("Angle bond (%i) at particle %i angle %lf partners: %i %i", bb.bid, bb.pid, angle, pid1, pid2);

                    printf(" pair bonds:");
                    // Search for the opposite and neighboring pair bonds
                    std::array<std::pair<std::optional<BondStore>, double>, 3> obs = {{
                        { find_pair_bond(fbs, pid1, pid2), pdist(s, pid1, pid2) },
                        { find_pair_bond(fbs, pid, pid1), pdist(s, pid, pid1) },
                        { find_pair_bond(fbs, pid, pid2), pdist(s, pid, pid2) },
                    }};

                    for (const auto [ob, dist]: obs) {
                        if (ob) {
                            if (dist >= 2.0)
                                printf(" broken (%i: %lf)", ob->bid, dist);
                            else
                                printf(" valid (%i: %lf)", ob->bid, dist);
                        } else {
                                printf(" none ()");
                        }
                    }

                    printf(" angle bonds:");
                    std::array<std::array<particle_id, 3>, 2> triangles = {{
                        {{ pid1, pid, pid2 }},
                        {{ pid2, pid, pid1 }},
                    }};

                    for (const auto [t1, t2, t3]: triangles) {
                        auto ob = find_angle_bond(fbs, t1, t2, t3);
                        if (ob) {
                            auto angle = calc_bond_angle(s, *ob);
                            if (is_broken_angle(angle))
                                printf(" broken (%i: %lf)", ob->bid, angle);
                            else
                                printf(" valid (%i: %lf)", ob->bid, angle);
                        } else {
                                printf(" none ()");
                        }
                    }

                    putchar('\n');
                }
            }
        }
    }
}

void eusage(const char *argv0)
{
    std::fprintf(stderr, "Usage: %s SNAP-PREFIX FUNC\n", argv0);
    std::exit(1);
}

int main(int argc, char **argv)
{
    if (argc < 2)
        eusage(*argv);
    
    auto s = snapshot{argv[1]};
    auto bs = BondingStructure{s.npart()};
    auto fbs = FullBondStorage{};

    bool store_full = false;
    if (argc > 2 && std::string(argv[2]) == "--full-storage") {
        fbs.resize(s.npart());
        store_full = true;
        argv++; argc--;
    }

    auto particle_callback = [](auto){};
    auto bond_callback = [&bs, &s, store_full, &fbs](Bond b){
        for (int i = 0; i < b.npartners; ++i)
            bs.add_bond(b.pid, b.partner_ids[i], pdist(s, b.pid, b.partner_ids[i]));
        if (b.npartners == 1) {
            bs.add_hbond(b.pid, b.partner_ids[0], pdist(s, b.pid, b.partner_ids[0]));
        }
        if (b.npartners == 2) {
            double angle = calc_bond_angle(s, b);
            bs.add_angle(b.pid, b.partner_ids[0], b.partner_ids[1], angle);
#ifdef REPORT_BROKEN_ANGLES
            if (is_broken_angle(angle))
                printf("Broken angle: %lf %i\n", angle, b.bid);
#endif
        }

        if (store_full)
            fbs[b.pid].push_back(b);
    };

    snapshot_iter(s, particle_callback, bond_callback);

    auto aggs = bs.agglomerates();


    for (char **argi = argv + 2; argi < argv + argc; ++argi) {
        auto arg = std::string{*argi};

        if (arg == "--sizes") {
            print_agglomerate_sizes(aggs);
        } else if (arg == "--size-with-maxdist") {
            print_agglomerate_sizes_with_maxdists(s, aggs);
        } else if (arg == "--size-with-maxbl") {
            print_agglomerate_sizes_with_maxbl(bs, aggs);
        } else if (arg == "--size-with-maxangle") {
            print_agglomerate_sizes_with_maxangle(bs, aggs);
        } else if (arg == "--size-with-all") {
            print_agglomerate_sizes_with_all(s, bs, aggs);
        } else if (arg == "--size-hist") {
            print_agglomerate_sizes_histogram(aggs);
        } else if (arg == "--print") {
            auto sz = strtoul(*(argi + 1), nullptr, 10);
            print_agglomerate_of_size(s, aggs, sz, true, true, true);
            argi++;
        } else if (arg == "--print-to-files") {
            std::vector<size_t> wanted;
            for (argi++; *argi && (*argi)[0] != '-'; argi++) {
                auto sz = strtoul(*argi, nullptr, 10);
                wanted.push_back(sz);
            }
            print_agglomerates_of_sizes(s, bs, aggs, wanted);
        } else if (arg == "--all-to-files") {
            print_agglomerates_to_files(s, bs, aggs);
        } else if (arg == "--print-pos") {
            auto sz = strtoul(*(argi + 1), nullptr, 10);
            print_agglomerate_of_size(s, aggs, sz, true, false);
            argi++;
        } else if (arg == "--print-vel") {
            auto sz = strtoul(*(argi + 1), nullptr, 10);
            print_agglomerate_of_size(s, aggs, sz, false, true);
            argi++;
        } else if (arg == "--print-raw") {
            ensure(store_full);
            auto sz = strtoul(*(argi + 1), nullptr, 10);
            print_agglomerate_raw(s, fbs, aggs, sz);
            argi++;
        } else if (arg == "--print-broken-bondlength") {
            auto sz = strtoul(*(argi + 1), nullptr, 10);
            print_agglomerate_of_size_broken_bonds(bs, s, aggs, sz);
            argi++;
        } else if (arg == "--print-broken-angle") {
            auto sz = strtoul(*(argi + 1), nullptr, 10);
            print_agglomerate_of_size_broken_angle_bonds(bs, s, aggs, sz);
            argi++;
        } else if (arg == "--print-all-broken-bonds") {
            print_all_broken_bonds(bs, s, aggs);
        } else if (arg == "--print-all-broken-bonds-verbose") {
            print_all_broken_bonds_verbose(bs, s, aggs);
        } else if (arg == "--df") {
            print_agglomerate_dfs(s, aggs);
        } else if (arg == "--check-angle-bonds") {
            ensure(store_full);
            check_angle_bonds(bs, s, aggs, fbs);
        } else if (arg == "--track-particle") {
            ensure(store_full);
            auto pid = atoi(*(argi + 1));
            track_particle(bs, s, aggs, fbs, pid);
            argi++;
        } else if (arg == "--print-agglomerate-of-pid") {
          auto pid = atoi(*(argi + 1));
          print_agglomerate_of_pid(s, aggs, pid);
          argi++;
        } else if (arg == "--print-pos-broken-by-angle") {
            print_broken_angle_bonds_with_positions(bs, s, aggs);
        } else {
            std::fprintf(stderr, "Unrecognized option: %s\n", *argi);
            std::exit(1);
        }
    };


}
