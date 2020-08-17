// For license details see LICENSE.

#include <cmath>
#include <cstdio>
#include <getopt.h>
#include <map>
#include <optional>
#include <string>

#include "bonding_structure.hpp"
#include "cfile.hpp"
#include "df.hpp"
#include "snapshot.hpp"
#include "uf.hpp"

std::vector<Vec3d> ids_to_poss_copy(const snapshot &s, const Agglomerate &agg) {
    std::vector<Vec3d> poss;
    poss.reserve(agg.size());
    for (const particle_id pid : agg) {
        const span3d p = s.pos_of_part(pid);
        poss.emplace_back(Vec3d{p[0], p[1], p[2]});
    }
    return poss;
}

void print_agglomerates_to_files(const snapshot &s,
                                 const std::vector<Agglomerate> &aggs) {
    std::map<int, int> num;
    for (const auto &agg : aggs) {
        auto [ii, unseen] = num.emplace(agg.size(), 0);
        ii->second++;
        auto fn =
            std::to_string(agg.size()) + "_POS_" + std::to_string(ii->second);
        if (auto f = cfile(fn, "w")) {
            for (const particle_id pid : agg) {
                const auto pos = s.pos_of_part(pid);
                std::fprintf(f, "%lf %lf %lf\n", pos[0], pos[1], pos[2]);
            }
        }
    }
}

void print_agglomerate_dfs(const snapshot &s,
                           const std::vector<Agglomerate> &aggs, double box_l,
                           double sigma) {
    auto calc_df_radog = [&s, sigma, box_l](const Agglomerate &agg) {
        return calc_df(ids_to_poss_copy(s, agg), box_l, sigma);
    };

    for (const auto &agg : aggs) {
        if (agg.size() < 15)
            continue;
        auto [radog, df_radog] = calc_df_radog(agg);
        std::printf("%zu %lf %lf\n", agg.size(), radog, df_radog);
    }
}

void eusage(const char *argv0) {
    std::fprintf(stderr, "Usage: %s [OPTIONS...] SNAP-PREFIX\n", argv0);

    constexpr const char *h =
        "Options:\n"
        " --sigma SIGMA | -s SIGMA   Set sigma for Df calculation\n"
        " --box_l BOX | -b BOX       Set box size for Df calculation\n"
        "Modi:\n"
        " --print-all-to-files | -p  Print all agglomerates to files named\n"
        "                              <NPART>_POS_<IDX>\n"
        " --df | -d                  Calculate fractal dimensions\n"
        "";
    std::fputs(h, stderr);
    std::exit(1);
}

static const option longopts[] = {
    {"sigma", required_argument, 0, 's'},
    {"box_l", required_argument, 0, 'b'},
    {"print-all-to-files", no_argument, 0, 'p'},
    {"df", no_argument, 0, 'd'},
    {0, 0, 0, 0},
};

struct {
    bool dflag = false;
    bool pflag = false;
    std::optional<double> sigma;
    std::optional<double> box_l;
} cmd_params;

void parse_args(int argc, char **argv) {
    int optidx;
    for (int c = 0; c != -1;
         c = getopt_long(argc, argv, "", longopts, &optidx)) {
        switch (c) {
        case 's':
            if (optarg)
                cmd_params.sigma = atof(optarg);
            break;
        case 'b':
            if (optarg)
                cmd_params.box_l = atof(optarg);
            break;
        case 'p':
            cmd_params.pflag = true;
            break;
        case 'd':
            cmd_params.dflag = true;
            break;
        case '?':
            eusage(*argv);
            break;
        default:
            break;
        }
    }

    /* Parameter sanity check */
    if (!cmd_params.pflag && !cmd_params.dflag) {
        std::fprintf(stderr, "No mode specified.\n");
        eusage(*argv);
    }

    if (cmd_params.dflag && (!cmd_params.box_l || !cmd_params.sigma)) {
        std::fprintf(stderr, "Df calculation requires --sigma and --box_l.\n");
        eusage(*argv);
    }

    if (optind != argc - 1) {
        std::fprintf(stderr, "No snapshot specified.\n");
        eusage(*argv);
    }
}

int main(int argc, char **argv) {
    parse_args(argc, argv);

    auto s = snapshot{argv[optind]};
    auto bs = BondingStructure{s.npart()};

    auto particle_callback = [](auto) {};
    auto bond_callback = [&bs, &s](BondReference b) {
        for (int i = 0; i < b.npartners; ++i)
            bs.add_bond(b.pid, b.partner_ids[i]);
    };

    snapshot_iter(s, particle_callback, bond_callback);

    auto aggs = bs.agglomerates();

    if (cmd_params.pflag)
        print_agglomerates_to_files(s, aggs);

    if (cmd_params.dflag)
        print_agglomerate_dfs(s, aggs, *cmd_params.box_l, *cmd_params.sigma);
}