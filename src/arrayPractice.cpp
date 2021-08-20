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
#include "arrayPractice.h"

#define ELEMENTS 2817

// Function declerations
std::vector<double> readFile(std::string file_name);
template<std::size_t SIZE>
std::array<double, SIZE> deriv(std::array<double, SIZE>& data);
template<std::size_t SIZE>
std::array<double, SIZE> fastSmooth(std::array<double, SIZE>& yData, unsigned int smoothWidth);
template<std::size_t SIZE>
std::vector<double> val2ind(std::array<double, SIZE>& xData, double val);
template<std::size_t SIZE>
std::array<double, SIZE> forEach(const std::array<double, SIZE>& values, double modifier, double(*func)(double, double));

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
    std::array<double, ELEMENTS> yDeriv = deriv(yData);
    if (smoothWidth > 1) {
        yDeriv = fastSmooth(yDeriv, smoothWidth);
    }

    // find peaks
    double peakX, peakY;
    int groupIndex;
    std::vector<double> pIndex;
    std::vector<double> allHeights;
    std::vector<double> allLocations;
    std::vector<int> allIndexes;
    for (int j = 2 * round(smoothWidth/2) - 2; j < ELEMENTS - smoothWidth - 1; j++)
    {
        if (yDeriv[j] > 0 && yDeriv[j + 1] < 0)
        {
            if (yDeriv[j]-yDeriv[j + 1] > slopeThreshold)
            {
                std::array<double, peakGroup> xIndexes;
                std::array<double, peakGroup> yIndexes;
                for (int k = 0; k < peakGroup; k++)
                {
                    groupIndex = j + k - round(peakGroup / 2 + 1) + 2;
                    if (groupIndex < 1) {groupIndex = 1; }
                    if (groupIndex > ELEMENTS) {groupIndex = ELEMENTS; }
                    xIndexes[k] = xData[groupIndex];
                    yIndexes[k] = yData[groupIndex];
                }
                
                if (peakGroup < 3) {peakY = *std::max_element(yIndexes.begin(),yIndexes.end()); }
                else {peakY = std::accumulate(yIndexes.begin(), yIndexes.end(), 0.0)/yIndexes.size(); }

                pIndex = val2ind(yIndexes, peakY);
                peakX = xIndexes[pIndex[0]];

                if (peakY > ampThreshold)
                {
                    allIndexes.push_back(j);
                    allHeights.push_back(peakY);
                    allLocations.push_back(peakX);
                    // printf("%i, height: %f, location: %f\n", j, peakY, peakX);
                }
            }
        }
    }

    std::vector<int> filteredIndexes;
    std::vector<double> filteredHeights;
    std::vector<double> filteredLocations;

    for (int i = 1; i < allHeights.size() - 1; i++)
    {
        if (allHeights[i] - allHeights[i-1] > heightThreshold)
        {
            filteredIndexes.push_back(allIndexes[i]);
            filteredHeights.push_back(allHeights[i]);
            filteredLocations.push_back(allLocations[i]);
        }
    }
    
    int group = 0;
    std::vector<int> groupsVect;
    for (int i = 0; i < filteredIndexes.size() - 1; i++)
    {
        groupsVect.push_back(group);
        if (filteredLocations[i + 1] - filteredLocations[i] > locationThreshold)
        {
            group = group + 1;
        }
    }
    groupsVect.push_back(group);

    std::vector<int> indexes;
    std::vector<double> heights;
    std::vector<double> locations;
    std::vector<int> tempIndexes;
    std::vector<double> tempHeights;
    std::vector<double> tempLocations;
    
    int groups = *std::max_element(groupsVect.begin(), groupsVect.end()) + 1;
    for (int i = 0; i < groups; i++)
    {
        tempIndexes.clear();
        tempHeights.clear();
        tempLocations.clear();
        for (int j = 0; j < groupsVect.size(); j++)
        {
            if (groupsVect.at(j) == i)
            {
                tempIndexes.push_back(filteredIndexes[j]);
                tempHeights.push_back(filteredHeights[j]);
                tempLocations.push_back(filteredLocations[j]);
            }
        }
        int index = std::find(tempHeights.begin(), tempHeights.end(), *std::max_element(tempHeights.begin(),tempHeights.end())) - tempHeights.begin();
        indexes.push_back(tempIndexes[index]);
        heights.push_back(tempHeights[index]);
        locations.push_back(tempLocations[index]);
    }

    for (int k = 0; k < indexes.size(); k++)
    {
        std::cout << indexes.at(k) << ", " << heights.at(k) << ", " << locations.at(k) << std::endl;
    }

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

template<std::size_t SIZE>
std::array<double, SIZE> forEach(const std::array<double, SIZE>& values, double modifier, double(*func)(double, double))
{
    std::array<double, SIZE> modifiedArray;
    for (int i = 0; i < modifiedArray.size(); i++)
        modifiedArray[i] = func(values[i], modifier);

    return modifiedArray;
}

template<std::size_t SIZE>
std::vector<double> val2ind(std::array<double, SIZE>& yInd, double val)
{
    std::array<double, SIZE> sumArray = forEach(yInd, val, [](double a, double b) ->double {return a - b; });
    std::array<double, SIZE> diff = forEach(yInd, 0, [](double a, double b) ->double {return abs(a);});
    double minDiff = *std::min_element(diff.begin(), diff.end());
    diff = forEach(diff, minDiff, [](double a, double b) ->double {return a - b; });

    std::vector<double> indexes;
    for (int i = 0; i < diff.size(); i++)
    {
        if (diff[i] == 0) {indexes.push_back(i); }
    }

    return indexes;
}
