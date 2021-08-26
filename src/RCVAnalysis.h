#ifndef RCVANALYSIS_H
#define RCVANALYSIS_H

#define ELEMENTS 2817

/* Include Files */

/* Function Declarations */
extern "C" int main();
extern "C" int findPeaksWrapper(double* xData, 
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
extern "C" int findFeaturesWrapper(double* xData, 
    double* yData, 
    int* indexes, 
    const unsigned int smoothWidth, 
    double peakArea[2], 
    double fullHeight[2], 
    double FWHM[2], 
    double xValsFullHeight[4], 
    double xValsHalfHeight[4] );
extern "C" void readFileWrapper(char* file_name, double data[2817]);

#endif
