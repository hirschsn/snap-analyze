// For license details see LICENSE.

#include <map>
#include <iterator>
#include <iostream>
#include <cmath>
#include "box.hpp"
#include "stat.hpp"
#include "snapshot.hpp"

void eusage(const char *argv0)
{
    std::fprintf(stderr, "Usage: %s SNAP-PREFIX FUNC\n", argv0);
    std::exit(1);
}

enum Mode {
    Print,
    Bondlength,
    BondlengthHistogram,
    PrintBondlengths,
    Ekin,
    BondTypeHist,
    NModes,
};

static constexpr int nmodes() { return Mode::NModes; }

#include <map>
std::map<std::string, Mode> argmap = {
    {"--print", Mode::Print},
    {"--bondlen", Mode::Bondlength},
    {"--bondlen-hist", Mode::BondlengthHistogram},
    {"--print-bondlens", Mode::PrintBondlengths},
    {"--ekin", Mode::Ekin},
    {"--bond-type-hist", Mode::BondTypeHist},
};

std::array<bool, nmodes()> parse_argv(int argc, char **argv)
{
    
    std::array<bool, nmodes()> m;
    std::fill(std::begin(m), std::end(m), false);
    std::for_each(argv, argv + argc, [&m](char *arg){
        m[argmap[arg]] = true;
    });
    return m;
}

void print_histogram(const Histogram<double>& hist)
{
    std::cout << "Bond histogram:" << std::endl;
    static const double hist_incr = 0.1;
    double oldd = 0.5;
    for (double d = oldd + hist_incr; d <= hist.max() + hist_incr; d += hist_incr) {
        size_t np = std::count_if(hist.begin(), hist.end(), [oldd,d](double bl){ return bl < d && bl >= oldd; });
        auto perc = static_cast<double>(np) / hist.nsamples();
        std::cout << "[" << oldd << "," << d << "): " << np << "(" << perc * 100 << " %)" << std::endl;

        oldd = d;
    }
}

void print_bond_type_hist(const std::map<int, size_t>& hist)
{
    for (const auto& el: hist)
        std::cout << el.first << " " << el.second << std::endl;
}

int main(int argc, char **argv)
{
    if (argc < 2)
        eusage(*argv);

    auto s = snapshot{argv[1]};
    auto flags = parse_argv(argc - 2, argv + 2);

    if (flags[Mode::Print])
        std::cout << *argv << std::endl
                  << "Npart: " << s.npart() << std::endl
                  << " nproc: " << s.nproc() << std::endl;

    Histogram<double> bond_len;
    ScalarStatistics<double> ekin;
    std::map<int, size_t> bond_type_hist;

    auto particle_callback = [&flags, &s, &ekin](Particle p){
        if (flags[Mode::Print])
            std::printf("Particle %i (pos %lf %lf %lf), bond info:",
                        p.id, p.pos[0], p.pos[1], p.pos[2]);
        if (flags[Mode::Ekin])
            ekin.sample(kinetic_energy(s, p.id));
    };

    auto bond_callback = [&flags, &s, &bond_len, &bond_type_hist](Bond b){
        if (flags[Mode::Print]) {
            printf("<Bond npartners=%i particles=",  b.npartners);
            for (int i = 0; i < b.npartners; ++i) {
                printf(" %i", b.partner_ids[i]);
            }
        }

        if (b.npartners == 1) {
            auto dist = pdist(s, b.pid, b.partner_ids[0]);
            
            if (flags[Mode::Print])
                printf(" dist=%lf", dist);

            if (flags[Mode::Bondlength] || flags[Mode::BondlengthHistogram] || flags[Mode::PrintBondlengths])
                bond_len.sample(dist);
        }

        if (flags[Mode::Print])
            printf(">");


        if (flags[Mode::BondTypeHist]) {
            auto [it, unseen] = bond_type_hist.emplace(b.bid, static_cast<size_t>(0));
            it->second++;
        }
    };

    snapshot_iter(s, particle_callback, bond_callback);

    if (flags[Mode::Bondlength] || flags[Mode::Ekin])
        std::cout << argv[1];
    if (flags[Mode::Bondlength])
        std::cout << " BL " << bond_len;
    if (flags[Mode::Ekin]) {
        p_assert(ekin.nsamples() == s.npart());
        std::cout << " Ekin " << ekin;
    }
    if (flags[Mode::Bondlength] || flags[Mode::Ekin])
        std::cout << std::endl;
    
    if (flags[Mode::PrintBondlengths])
        std::copy(std::begin(bond_len), std::end(bond_len), std::ostream_iterator<double>(std::cout, "\n"));

    if (flags[Mode::BondlengthHistogram])
        print_histogram(bond_len);
    
    if (flags[Mode::BondTypeHist])
        print_bond_type_hist(bond_type_hist);
}