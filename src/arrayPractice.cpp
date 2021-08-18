#include <iostream>
#include <array>
#include <vector>
#include "stdio.h"
#include <stdlib.h>
#include <fstream>
#include <string.h>
#include <math.h>
#include <numeric>
#include "arrayPractice.h"

#define ELEMENTS 2817

// Function declerations
std::vector<double> readFile(std::string file_name);
template<std::size_t SIZE>
std::array<double, SIZE> deriv(std::array<double, SIZE>& data);
template<std::size_t SIZE>
std::array<double, SIZE> fastSmooth(std::array<double, SIZE>& yData, unsigned int smoothWidth);

int main()
{
    std::vector<double> xDataVect = readFile("C:/Users/menwst/Documents/CPP/arrayPractice/data/xData");
    std::vector<double> yDataVect = readFile("C:/Users/menwst/Documents/CPP/arrayPractice/data/yData");

    std::array<double, ELEMENTS> xData;
    std::copy(xDataVect.begin(), xDataVect.begin() + ELEMENTS, xData.begin());
    std::array<double, ELEMENTS> yData;
    std::copy(yDataVect.begin(), yDataVect.begin() + ELEMENTS, yData.begin());

    const double slopeThreshold=-1;
    const double ampThreshold=0;
    const unsigned int smoothWidth=5;
    const int peakGroup=3;
    const double heightThreshold=1;
    const double locationThreshold=0.02;
    
    // smooth data and find derivative
    std::array<double, ELEMENTS> deriv_yData = deriv(yData);
    std::array<double, ELEMENTS> smooth_yData = fastSmooth(deriv_yData, smoothWidth);

    // for (int i = 0; i < smooth_yData.size(); i++)
    //     printf("%f\n",smooth_yData[i]);

    std::cin.get();
    return 0;
}


std::vector<double> readFile(std::string file_name)
{
    std::vector<double> record;

    std::ifstream file;
    file.open(file_name);

    std::string column_one;

    while(getline(file, column_one))
    {
        record.push_back(atof(column_one.c_str()));
    }

    return record;
}


template<std::size_t SIZE>
std::array<double, SIZE> deriv(std::array<double, SIZE>& data)
{
    int numElements = data.size();
    std::array<double, SIZE> deriv;
    deriv[0] = data[1] - data[0];
    for(int j = 1; j < numElements - 1; j++)
    {
        deriv[j] = (data[j + 1] - data[j - 1]) / 2;
    }
    deriv[numElements - 1] = data[numElements - 1] - data[numElements - 2];

    return deriv;
}


template<std::size_t SIZE>
std::array<double, SIZE> fastSmooth(std::array<double, SIZE>& yData, unsigned int smoothWidth)
{
    std::array<double, SIZE> smooth_yData;
    std::array<double, SIZE> temp;
    double sumPoints = std::accumulate(yData.begin(), yData.begin() + smoothWidth, 0.0);
    int halfWidth = round(smoothWidth / 2);
    int numPoints = yData.size();
    for (int i = 0; i < numPoints - smoothWidth; i++)
    {
        temp[i + halfWidth] = sumPoints;
        sumPoints = sumPoints - yData[i];
        sumPoints = sumPoints + yData[i + smoothWidth];
    }
    // temp[numPoints - smoothWidth + halfWidth] = std::accumulate(yData.begin() + numPoints - smoothWidth, yData.begin() + numPoints - 1, 0.0);

    for (int i = 0; i < numPoints; i++)
    {
        smooth_yData[i] = temp[i] / smoothWidth;
    }

    return smooth_yData;
}
