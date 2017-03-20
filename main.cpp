#include <iostream>
#include "Miller.h"
using namespace std;

int main()
{
    // use cols > rows!


    MurtyMiller<double>::WeightMatrix c_ij(20, 30);
    typedef MurtyMiller<double>::Edges Edges;

    c_ij = MurtyMiller<double>::WeightMatrix::Random(20, 30);

    for ( size_t r = 0; r < c_ij.rows(); ++r )
        for ( size_t c = 0; c < c_ij.cols(); ++c )
            if ( c_ij(r, c) < 0 ) c_ij(r, c) = -c_ij(r, c);

    c_ij /= c_ij.maxCoeff();

    std::cerr << "c_ij = \n" << c_ij << std::endl;

    std::vector<Edges> solutions = MurtyMiller<double>::getMBestAssignments(c_ij, 100);

    for ( const auto & s : solutions )
    {
        for ( const auto & e : s )
            std::cerr << "(" << e.x << ", " << e.y << ") ";
        std::cerr << "sum = " << MurtyMiller<double>::objectiveFunctionValue(s) << std::endl;
    }

    return 0;
}

