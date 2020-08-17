// For license details see LICENSE.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <iterator>
#include <numeric>
#include <vector>

typedef std::array<double, 3> Vec3d;

Vec3d &operator+=(Vec3d &a, const Vec3d &b) {
    for (int i = 0; i < 3; ++i)
        a[i] += b[i];
    return a;
}

Vec3d &operator/=(Vec3d &x, double d) {
    for (int i = 0; i < 3; ++i)
        x[i] /= d;
    return x;
}

/** Squared distance between two vectors.
 */
double dist2(const Vec3d &a, const Vec3d &b) {
    double d2 = 0.0;
    for (int i = 0; i < 3; ++i) {
        double d = b[i] - a[i];
        d2 += d * d;
    }
    return d2;
}

/** Normalizes all positions s.t. pos[0] is at HALF_BOX_L in each of the
 * three spatial dimensions.
 */
void normalize_against_fist_pos(std::vector<Vec3d> &pos, double box_l) {
    const auto half_box_l = box_l / 2.;
    auto base = pos[0]; // Copy
    for (auto &p : pos) {
        for (int i = 0; i < 3; ++i) {
            p[i] = std::fmod(p[i] - (base[i] - half_box_l),
                             box_l); //  Box_l and box_l/2!
            if (p[i] < 0.0)
                p[i] += box_l;
        }
    }
}

/** Calculates the center of mass of a set of positions.
 */
Vec3d center_of_mass(const std::vector<Vec3d> &pos) {
    Vec3d com = {{0.0, 0.0, 0.0}};
    for (const auto &p : pos)
        com += p;
    com /= pos.size();
    return com;
}

/** Returns the slope of a linear regression.
 * [first1, last1) is the range of xs.
 * [first2, first2 + (last1 - first1)) is the range of ys.
 */
template <typename T, typename It1, typename It2>
T linregress(It1 first1, It1 last1, It2 first2) {
    size_t n = 0;
    T s_x = T(0), s_y = T(0), s_xx = T(0), s_xy = T(0);
    while (first1 != last1) {
        T xi = *first1++;
        T yi = *first2++;
        s_x += xi;
        s_y += yi;
        s_xx += xi * xi;
        s_xy += xi * yi;
        n++;
    }
    return (n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);
}

/** Calculate the Radius of gyration and the fractal dimension of a
 * single agglomerate.
 * Radius of gyration if the standard deviation of the particle positions
 * from the center of mass.
 * Fractal dimension is the power-law relationship of radius of gyration
 * to number of particles.
 */
std::pair<double, double>
calc_df(std::vector<Vec3d> pos /* take a deep copy, we modify it */,
        double box_l, double sigma) {
    normalize_against_fist_pos(pos, box_l);
    auto com = center_of_mass(pos);

    std::vector<double> dists(pos.size(), 0.0);
    for (size_t i = 0; i < pos.size(); ++i)
        dists[i] = std::sqrt(dist2(com, pos[i]));

    std::sort(std::begin(dists), std::end(dists));

    int min_dist = dists.front() / sigma;
    int max_dist = dists.back() / sigma;
    // Because of distance to center of mass, dists might not start with "0"
    size_t nsamples = max_dist - min_dist + 1;

    std::vector<size_t> ks(nsamples, 0);
    std::vector<double> dks(nsamples);
    std::vector<double> radogs(nsamples, 0.0), lradogs(nsamples, 0.0);

    for (auto d : dists) {
        int bin = std::floor(d / sigma - min_dist);
        ks[bin]++;
    }

    radogs[0] = std::sqrt(std::inner_product(std::begin(dists),
                                             std::begin(dists) + ks[0],
                                             std::begin(dists), 0.0) /
                          ks[0]);
    dks[0] = std::log(static_cast<double>(ks[0]));
    for (size_t i = 1; i < nsamples; ++i) {
        ks[i] = ks[i - 1] + ks[i];
        dks[i] = std::log(static_cast<double>(ks[i]));
        radogs[i] = std::sqrt(std::inner_product(std::begin(dists),
                                                 std::begin(dists) + ks[i],
                                                 std::begin(dists), 0.0) /
                              ks[i]);
    }

    std::transform(std::begin(radogs), std::end(radogs), std::begin(lradogs),
                   [](double d) { return std::log(d); });

    auto df_radog = linregress<double>(std::begin(lradogs), std::end(lradogs),
                                       std::begin(dks));

    return std::make_pair(radogs.back(), df_radog);
}
