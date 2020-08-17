// For license details see LICENSE.
 
#include <cmath>
#include <map>
#include <string>
#include <cstdio>

#include "snapshot.hpp"
#include "cfile.hpp"
#include "bonding_structure.hpp"
#include "df.hpp"
#include "uf.hpp"

template <typename Ret>
std::vector<Ret> __ids_to_poss_impl(const snapshot &s, const std::vector<int> &agg)
{
    std::vector<Ret> poss;
    poss.reserve(agg.size());
    std::transform(agg.begin(), agg.end(), std::back_inserter(poss), [&](particle_id i){ return s.pos_of_part(i);});
    return poss;
}

std::vector<span3d> ids_to_poss(const snapshot &s, const std::vector<particle_id> &agg)
{
    return __ids_to_poss_impl<span3d>(s, agg);
}

std::vector<Vec3d> ids_to_poss_copy(const snapshot &s, const std::vector<int> &agg)
{
    return __ids_to_poss_impl<Vec3d>(s, agg);
}

void print_agglomerates_to_files(const snapshot& s, const std::vector<std::vector<int>>& aggs)
{
    std::map<int, int> num;
    for (const auto& agg: aggs) {
        auto [ii, unseen] = num.emplace(agg.size(), 0);
        ii->second++;
        auto pos = ids_to_poss(s, agg);
        auto fn = std::to_string(agg.size()) + "_POS_" + std::to_string(ii->second);
        if (auto f = cfile(fn.c_str(), "w")) {
            for (size_t i = 0; i < agg.size(); ++i) {
                std::fprintf(f, "%lf %lf %lf\n", pos[i][0], pos[i][1], pos[i][2]);
            }
        }
    }
}

void print_agglomerate_dfs(const snapshot& s, const std::vector<std::vector<int>>& aggs)
{
    auto calc_df_radog = [&s](const std::vector<int> &agg){
        return calc_df(ids_to_poss_copy(s, agg), s.box_l);
    };

    for (const auto &agg: aggs) {
        if (agg.size() < 15)
            continue;
        auto [radog, df_radog] = calc_df_radog(agg);
        std::printf("%zu %lf %lf\n", agg.size(), radog, df_radog);
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
    
    // TODO: Pass via argument
    const double box_l = 800.;
    auto s = snapshot{box_l, argv[1]};
    auto bs = BondingStructure{s.npart()};

    auto particle_callback = [](auto){};
    auto bond_callback = [&bs, &s](BondReference b){
        for (int i = 0; i < b.npartners; ++i)
            bs.add_bond(b.pid, b.partner_ids[i]);
    };

    snapshot_iter(s, particle_callback, bond_callback);

    auto aggs = bs.agglomerates();


    for (char **argi = argv + 2; argi < argv + argc; ++argi) {
        auto arg = std::string{*argi};

        if (arg == "--all-to-files") {
            print_agglomerates_to_files(s, aggs);
        } else if (arg == "--df") {
            print_agglomerate_dfs(s, aggs);
        } else {
            std::fprintf(stderr, "Unrecognized option: %s\n", *argi);
            std::exit(1);
        }
    };


}
