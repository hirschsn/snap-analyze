// For license details see LICENSE.
//
// This program calculates the radius of gyration and the fractal dimension
// of a set of 3d particle positions from standard input.
//

#include <iostream>
#include "df.hpp"

int main(void)
{
    // TODO: Pass via argument
    const double box_l = 800.;
    
    std::vector<Vec3d> vs;
    while (std::cin) {
        Vec3d v;
        std::cin >> v[0] >> v[1] >> v[2];
        vs.push_back(std::move(v));
    }

    auto [radog, df_radog] = calc_df(vs, box_l);
    std::cout << "Radog Df_radog" << std::endl;
    std::cout << radog << " " << df_radog << std::endl;
}