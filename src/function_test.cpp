#include <iostream>
#include <array>
#include <vector>
#include "stdio.h"
#include <stdlib.h>
#include <fstream>
#include <string.h>
#include <math.h>
#include <numeric>

int main()
{
    // sum array elements with std::accumulate
    std::array<double, 4>yData = {1.0, 2.0, 1.5, 2.0};
    // double yData[4] = {1.0, 2.0, 1.5, 2.0};
    double sumPoints = std::accumulate(yData.begin() + 1, yData.end(), 0.0);
    printf("%f\n", sumPoints);

    // sum array elements with for loop
    double sum_of_elems;
    for(int it = 0; it < 3; ++it)
    {
        sum_of_elems += yData[it];
    }
    printf("%f\n", sum_of_elems);

    double halfWidth = round(2.5);
    printf("%f\n", halfWidth);

    for(auto val : yData)
    {
         printf("%f, ", val);
    }

    std::cin.get();
    return 0;
}