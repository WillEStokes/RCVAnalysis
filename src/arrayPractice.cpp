#include <iostream>
#include <array>
#include <vector>
#include "stdio.h"
#include <stdlib.h>
#include <fstream>
#include <string.h>
#include <math.h>
#include "arrayPractice.h"

#define ELEMENTS 2817

// Function declerations
std::vector<double> readFile(std::string file_name);
template<std::size_t SIZE>
std::array<double, SIZE> deriv(std::array<double, SIZE>& data);

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
    const int smoothWidth=5;
    const int peakGroup=3;
    const double heightThreshold=1;
    const double locationThreshold=0.02;
    
    // smooth data and find derivative
    std::array<double, ELEMENTS> deriv_yData = deriv(yData);

    for (int i=0; i < deriv_yData.size(); i++)
    {
        printf("%f\n",deriv_yData[i]);
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
        deriv[j] = (data[j+1]-data[j-1]) / 2;
    }
    deriv[numElements] = data[numElements] - data[numElements-1];

    return deriv;
}


// function SmoothY=fastsmooth(Y,smoothwidth)
// SumPoints=sum(Y(1:smoothwidth));
// s=zeros(size(Y));
// halfw=round(smoothwidth/2);
// L=length(Y);
// for k=1:L-smoothwidth
//     s(k+halfw-1)=SumPoints;
//     SumPoints=SumPoints-Y(k);
//     SumPoints=SumPoints+Y(k+smoothwidth);
// end
// s(k+halfw)=sum(Y(L-smoothwidth+1:L));
// SmoothY=s./smoothwidth;
// end