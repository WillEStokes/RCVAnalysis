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
#include <chrono> // for timing
#include <unistd.h> // for sleep
#include <ctime>
#include "RCVAnalysis.h"

// Function declerations
void findPeaksWrapper(double* xData, 
    double* yData, 
    const double slopeThreshold, 
    const double ampThreshold, 
    const unsigned int smoothWidth, 
    const unsigned int peakGroup, 
    const double heightThreshold, 
    const double locationThreshold, 
    int indexes[2], 
    double heights[2], 
    double locations[2] );
void findFeaturesWrapper(double* xData, 
    double* yData, 
    int* indexes, 
    const unsigned int smoothWidth, 
    double peakArea[2], 
    double fullHeight[2], 
    double FWHM[2], 
    double xValsFullHeight[4], 
    double xValsHalfHeight[4] );
void findPeaks(std::array<double, ELEMENTS> xData, 
    std::array<double, ELEMENTS> yData, 
    const double slopeThreshold, 
    const double ampThreshold, 
    const unsigned int smoothWidth, 
    const unsigned int peakGroup, 
    const double heightThreshold, 
    const double locationThreshold, 
    std::array<int, 2>& indexes, 
    std::array<double, 2>& heights, 
    std::array<double, 2>& locations );
void findFeatures(std::array<double, ELEMENTS> xData, 
    std::array<double, ELEMENTS> yData, 
    std::array<int, 2> indexes, 
    const unsigned int smoothWidth, 
    std::array<double, 2>& peakArea, 
    std::array<double, 2>& fullHeight, 
    std::array<double, 2>& FWHM, 
    std::array<std::array<double, 2>, 2>& xValsFullHeight, 
    std::array<std::array<double, 2>, 2>& xValsHalfHeight );
void readFile(std::string file_name, std::vector<double>& data);
template<std::size_t SIZE>
std::array<double, SIZE> deriv(std::array<double, SIZE>& data);
template<std::size_t SIZE>
std::array<double, SIZE> fastSmooth(std::array<double, SIZE>& yData, unsigned int smoothWidth);
template<std::size_t SIZE>
std::vector<double> val2ind(std::array<double, SIZE>& xData, int numPeaks, double val);
template<std::size_t SIZE>
std::array<double, SIZE> forEach(const std::array<double, SIZE>& values, double modifier, double(*func)(double, double));
double trapz(std::array<double, ELEMENTS> xData, std::array<double, ELEMENTS> yData, int beginInd, int endInd);

int main()
{
    std::vector<double> xDataVect;
    std::vector<double> yDataVect;
    readFile("C:/Users/menwst/Documents/CPP/RCVAnalysis/data/xData", xDataVect);
    readFile("C:/Users/menwst/Documents/CPP/RCVAnalysis/data/yData", yDataVect);

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
    std::array<int, 2> indexes;
    std::array<double, 2> heights;
    std::array<double, 2> locations;
    
    // auto t1 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    
    findPeaks(xData, yData, slopeThreshold, ampThreshold, smoothWidth, peakGroup, heightThreshold, locationThreshold, indexes, heights, locations);
    
    // auto t2 = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    
    std::array<double, 2> peakArea;
    std::array<double, 2> fullHeight;
    std::array<double, 2> FWHM;
    std::array<std::array<double, 2>, 2> xValsFullHeight;
    std::array<std::array<double, 2>, 2> xValsHalfHeight;

    findFeatures(xData, yData, indexes, smoothWidth, peakArea, fullHeight, FWHM, xValsFullHeight, xValsHalfHeight);

    // auto t2 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    // std::cout << "Elapsed time: " << t2 - t1 << " ns" << "; t2 = " << t2 << " ns; t1 = " << t1 << " ns" << std::endl;

    for (int k = 0; k < indexes.size(); k++)
    {
        std::cout << indexes.at(k) << ", " << heights.at(k) << ", " << locations.at(k) << std::endl;
    }

    std::cout << "peak 1 area: " << peakArea[0] << "; peak 2 area: " << peakArea[1] << std::endl;
    
    std::cin.get();
    return 0;
}

void findPeaksWrapper(double* xData, 
    double* yData, 
    const double slopeThreshold, 
    const double ampThreshold, 
    const unsigned int smoothWidth, 
    const unsigned int peakGroup, 
    const double heightThreshold, 
    const double locationThreshold, 
    int indexes[2], 
    double heights[2], 
    double locations[2] )
{
    std::array<double, ELEMENTS> xDataArray;
    memcpy(xDataArray.data(), xData, ELEMENTS*sizeof(double));
    std::array<double, ELEMENTS> yDataArray;
    memcpy(yDataArray.data(), yData, ELEMENTS*sizeof(double));
    std::array<int, 2> indexesArray;
    memcpy(indexesArray.data(), indexes, 2*sizeof(int));
    std::array<double, 2> heightsArray;
    memcpy(heightsArray.data(), heights, 2*sizeof(double));
    std::array<double, 2> locationsArray;
    memcpy(locationsArray.data(), locations, 2*sizeof(double));

    findPeaks(xDataArray, yDataArray, slopeThreshold, ampThreshold, smoothWidth, peakGroup, heightThreshold, locationThreshold, indexesArray, heightsArray, locationsArray);

    for (int i = 0; i < 2; i++)
    {
        indexes[i] = indexesArray[i];
        heights[i] = heightsArray[i];
        locations[i] = locationsArray[i];
    }
}

void findFeaturesWrapper(double* xData, 
    double* yData, 
    int* indexes, 
    const unsigned int smoothWidth, 
    double peakArea[2], 
    double fullHeight[2], 
    double FWHM[2], 
    double xValsFullHeight[4], 
    double xValsHalfHeight[4] )
{
    std::array<double, ELEMENTS> xDataArray;
    memcpy(xDataArray.data(), xData, ELEMENTS*sizeof(double));
    std::array<double, ELEMENTS> yDataArray;
    memcpy(yDataArray.data(), yData, ELEMENTS*sizeof(double));
    std::array<int, 2> indexesArray;
    memcpy(indexesArray.data(), indexes, 2*sizeof(int));
    std::array<double, 2> peakAreaArray;
    memcpy(peakAreaArray.data(), peakArea, 2*sizeof(double));
    std::array<double, 2> fullHeightArray;
    memcpy(fullHeightArray.data(), fullHeight, 2*sizeof(double));
    std::array<double, 2> FWHMArray;
    memcpy(FWHMArray.data(), FWHM, 2*sizeof(double));
    std::array<std::array<double, 2>, 2> xValsFullHeightArray;
    memcpy(xValsFullHeightArray.data(), xValsFullHeight, 4*sizeof(double));
    std::array<std::array<double, 2>, 2> xValsHalfHeightArray;
    memcpy(xValsHalfHeightArray.data(), xValsHalfHeight, 4*sizeof(double));

    findFeatures(xDataArray, yDataArray, indexesArray, smoothWidth, peakAreaArray, fullHeightArray, FWHMArray, xValsFullHeightArray, xValsHalfHeightArray);

    int it = 0;
    for (int i = 0; i < 2; i++)
    {
        peakArea[i] = peakAreaArray[i];
        fullHeight[i] = fullHeightArray[i];
        FWHM[i] = FWHMArray[i];
        for (int j = 0; j < 2; j++)
        {
            xValsFullHeight[it] = xValsFullHeightArray[i][j];
            xValsHalfHeight[it] = xValsHalfHeightArray[i][j];
            it = it + 1;
        }
    }

}

void findPeaks(std::array<double, ELEMENTS> xData, 
    std::array<double, ELEMENTS> yData, 
    const double slopeThreshold, 
    const double ampThreshold, 
    const unsigned int smoothWidth, 
    const unsigned int peakGroup, 
    const double heightThreshold, 
    const double locationThreshold, 
    std::array<int, 2>& indexes, 
    std::array<double, 2>& heights, 
    std::array<double, 2>& locations )
{
    // smooth data and find derivative
    std::array<double, ELEMENTS> yDeriv = deriv(yData);
    if (smoothWidth > 1) {
        yDeriv = fastSmooth(yDeriv, smoothWidth);
    }

    // find peaks
    double peakX, peakY;
    int groupIndex;
    int numPeaks;
    std::array<double, ELEMENTS> xIndexes;
    std::array<double, ELEMENTS> yIndexes;
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
                numPeaks = 0;
                xIndexes.fill(0);
                yIndexes.fill(0);

                for (int k = 0; k < peakGroup; k++)
                {
                    groupIndex = j + k - round(peakGroup / 2 + 1) + 2;
                    if (groupIndex < 1) {groupIndex = 1; }
                    if (groupIndex > ELEMENTS) {groupIndex = ELEMENTS; }
                    xIndexes[k] = xData[groupIndex];
                    yIndexes[k] = yData[groupIndex];
                    numPeaks = numPeaks + 1;
                }
                
                if (peakGroup < 3) {peakY = *std::max_element(yIndexes.begin(), yIndexes.begin() + numPeaks); }
                else {peakY = std::accumulate(yIndexes.begin(), yIndexes.begin() + numPeaks, 0.0) / numPeaks; }

                pIndex = val2ind(yIndexes, numPeaks, peakY);
                peakX = xIndexes[pIndex[0]];

                if (peakY > ampThreshold)
                {
                    allIndexes.push_back(j);
                    allHeights.push_back(peakY);
                    allLocations.push_back(peakX);
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

    std::vector<int> tempIndexes;
    std::vector<double> tempHeights;
    std::vector<double> tempLocations;
    
    // int groups = *std::max_element(groupsVect.begin(), groupsVect.end()) + 1;
    // for (int i = 0; i < groups; i++)
    for (int i = 0; i < 2; i++)
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
        indexes[i] = tempIndexes[index];
        heights[i] = tempHeights[index];
        locations[i] = tempLocations[index];
    }

}

void findFeatures(std::array<double, ELEMENTS> xData, 
    std::array<double, ELEMENTS> yData, 
    std::array<int, 2> indexes, 
    const unsigned int smoothWidth, 
    std::array<double, 2>& peakArea, 
    std::array<double, 2>& fullHeight, 
    std::array<double, 2>& FWHM, 
    std::array<std::array<double, 2>, 2>& xValsFullHeight, 
    std::array<std::array<double, 2>, 2>& xValsHalfHeight )
{
    yData = fastSmooth(yData, smoothWidth);
    std::array<double, ELEMENTS> yDeriv = deriv(yData);

    int xIndsFullHeight[2][2];
    double yValsFullWidth[2][2];
    for (int j = 0; j < 2; j++ )
    {
        for (int i = indexes[j] - 1; i >= 0; i-- )
        {
            if (yDeriv[i] < 0)
            {
                xValsFullHeight[j][0] = xData[i];
                yValsFullWidth[j][0] = yData[i];
                xIndsFullHeight[j][0] = i;
                break;
            }
        }
        for (int i = indexes[j] + 1; i < ELEMENTS; i++ )
        {
            if (yDeriv[i] > 0)
            {
                xValsFullHeight[j][1] = xData[i];
                yValsFullWidth[j][1] = yData[i];
                xIndsFullHeight[j][1] = i;
                break;
            }
        }
    }

    for (int j = 0; j < 2; j++ )
    {
        fullHeight[j] = yData[indexes[j]] - (yValsFullWidth[j][0] + (yValsFullWidth[j][1] - yValsFullWidth[j][0]) / 2 );
        for (int i = indexes[j] - 1; i >= 0; i-- )
        {
            if (yData[i] < yData[indexes[j]] - fullHeight[j] / 2 )
            {
                xValsHalfHeight[j][0] = xData[i];
                break;
            }
        }
        for (int i = indexes[j] + 1; i < ELEMENTS; i++ )
        {
            if (yData[i] < yData[indexes[j]] - fullHeight[j] / 2 )
            {
                xValsHalfHeight[j][1] = xData[i];
                break;
            }
        }
        FWHM[j] = xValsHalfHeight[j][1] - xValsHalfHeight[j][0];
    }

    for (int j = 0; j < 2; j++ )
    {
        peakArea[j] = trapz(xData, yData, xIndsFullHeight[j][0], xIndsFullHeight[j][1]);
    }
}

void readFile(std::string file_name, std::vector<double>& data)
{
    std::ifstream file;
    file.open(file_name);

    std::string column_one;

    while(getline(file, column_one))
    {
        data.push_back(atof(column_one.c_str()));
    }
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
std::vector<double> val2ind(std::array<double, SIZE>& yInd, int numPeaks, double val)
{
    std::array<double, SIZE> sumArray = forEach(yInd, val, [](double a, double b) ->double {return a - b; });
    std::array<double, SIZE> diff = forEach(yInd, 0, [](double a, double b) ->double {return abs(a);});
    double minDiff = *std::min_element(diff.begin(), diff.begin() + numPeaks);
    diff = forEach(diff, minDiff, [](double a, double b) ->double {return a - b; });

    std::vector<double> indexes;
    for (int i = 0; i < numPeaks; i++)
    {
        if (diff[i] == 0) {indexes.push_back(i); }
    }

    return indexes;
}

double trapz(std::array<double, ELEMENTS> xData, std::array<double, ELEMENTS> yData, int beginInd, int endInd)
{
    // double dx = (xData[endInd] - xData[beginInd]) / (endInd - beginInd + 1);
    double dx = xData[1] - xData[0];
    
    double peakArea = dx * yData[beginInd] / 2;
    for (int i = beginInd + 1; i < endInd; i++)
    {
        peakArea = peakArea + dx * yData[i];
    }
    peakArea = peakArea + dx * yData[endInd] / 2;

    return peakArea;
}
