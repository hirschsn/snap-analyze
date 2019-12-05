
#include <iostream>
#include "df.hpp"

std::istream& operator>>(std::istream& is, Vec3d& x)
{
    is >> x[0] >> x[1] >> x[2];
    return is;
}

int main()
{
    //std::vector<Vec3d> pos = {
    //    {{1.0, 1.0, 1.0}},
    //    {{2.0, 2.0, 2.0}},
    //    {{2.0, 3.0, 3.0}},
    //    {{2.0, 4.0, 4.0}},
    //    {{2.0, 5.0, 5.0}},
    //    {{2.0, 6.0, 6.0}},
    //    {{2.0, 6.0, 3.0}},
    //    {{2.0, 5.0, 2.0}},
    //    {{2600.1, 2600.1, 2600.1}},
    //};

    std::vector<Vec3d> pos;

    while (std::cin) {
        Vec3d v;
        std::cin >> v;
        pos.push_back(v);
    }

    auto [radog, df_radog] = calc_df(std::move(pos));

    //std::cout << "Radog Df_radog" << std::endl;
    std::cout << radog << " " << df_radog << std::endl;
}