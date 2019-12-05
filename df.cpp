

#include <iostream>
#include "df.hpp"

int main(int argc, char **argv)
{
    if (argc != 2)
        exit(1);
    
    std::vector<Vec3d> vs;
    while (std::cin) {
        Vec3d v;
        std::cin >> v[0] >> v[1] >> v[2];
        vs.push_back(std::move(v));
    }

    auto [radog, df_radog] = calc_df(vs);
    std::cout << "Radog Df_radog" << std::endl;
    std::cout << radog << " " << df_radog << std::endl;
}