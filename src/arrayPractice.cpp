#include <iostream>
#include <array>
#include <vector>
#include "stdio.h"
#include <stdlib.h>
#include <fstream>
#include <string.h>
#include <math.h>
#include "arrayPractice.h"

std::vector<double> readFile(std::string file_name);

int main()
{
    std::vector<double> xData = readFile("C:/Users/menwst/Documents/CPP/arrayPractice/data/xData");
    std::vector<double> yData = readFile("C:/Users/menwst/Documents/CPP/arrayPractice/data/yData");
    
    // for (int i=0; i < xData.size(); i++)
    // {
        // printf("%f\n",xData[i]);
    // }

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
