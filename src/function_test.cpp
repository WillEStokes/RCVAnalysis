#include <iostream>
#include <array>
#include <vector>
#include "stdio.h"
#include <stdlib.h>
#include <fstream>
#include <string.h>
#include <math.h>
#include <numeric>
#include <algorithm>
#include <functional>
#include <chrono>

template<std::size_t SIZE>
std::array<double, SIZE> ForEach(const std::array<double, SIZE>& values, double modifier, double(*func)(double, double))
{
    std::array<double, SIZE> modifiedArray;
    for (int i = 0; i < modifiedArray.size(); i++)
        modifiedArray[i] = func(values[i], modifier);

    return modifiedArray;
}

int main()
{
    uint8_t newAddress = 104;
    // for (int i = 0; i < 2; i++) {
    // std::string tmp = std::to_string(newAddress);
    // char const *newAddressChar = tmp.c_str();
    // std::cout << tmp << std::endl;
    // std::cout << newAddressChar << std::endl;
    // std::cout << newAddressChar[0] << std::endl;
    // std::cout << newAddressChar[1] << std::endl;
    // std::cout << newAddressChar[2] << std::endl;
    // newAddress++;
    // }

    uint8_t value = static_cast<uint8_t>(newAddress);
    char dig[3] = { NULL, NULL, NULL};
    int x = 0;
    while (value > 0 && x < 3) {
        dig[x] = value % 10;
        value /= 10;
        x++;
    }
    // std::cout << dig[0] << dig[1] << dig[2] << std::endl;
    // printf("%i", dig[2]);
    // printf("%i", dig[1]);
    // printf("%i\n", dig[0]);
    // printf("x: %i\n", x);
    // printf("%i\n", value);

    if (dig[1] == NULL && dig[2] == NULL) {
        printf("%i\n", dig[0]);
    } else if (dig[2] == NULL) {
        printf("%i", dig[1]);
        printf("%i\n", dig[0]);
    } else {
        printf("%c", char(dig[2] + 48));
        printf("%i", dig[1]);
        printf("%i\n", dig[0]);
    }

    float time = 0.01;
    printf("chrono time: %i\n", std::chrono::microseconds(static_cast<long int>(time * 1000 / 2.0f)));
    
    const std::array<double, 4>yData = {1.0, 2.0, -1.5, 2.0};
    
    // index of element in vector
    int index = std::find(yData.begin(), yData.end(), -1.5) - yData.begin();
    printf("index: %i\n", index);

    // do maths to array elements with lambda function pointer
    std::array<double, 4> sumArray = ForEach(yData, 1.0, [](double a, double b) ->double {return a + b;});

    // abs of array
    std::array<double, 4> absArray = ForEach(yData, 1.0, [](double a, double b) ->double {return abs(a);});

    // min element of array
    double minElement = *std::min_element(yData.begin(), yData.end());
    printf("minELement: %f\n", minElement);

    for(int i = 0; i < 4; ++i)
        printf("%f, ", sumArray[i]);
    
    printf("\n");
    for(int i = 0; i < 4; ++i)
        printf("%f, ", absArray[i]);

    // sum array elements with std::accumulate
    double sumPoints = std::accumulate(yData.begin() + 1, yData.end(), 0.0);
    printf("\n%f\n", sumPoints);

    // sum array elements with for loop
    double sum_of_elems;
    for(int i = 0; i < 3; ++i)
    {
        sum_of_elems += yData[i];
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