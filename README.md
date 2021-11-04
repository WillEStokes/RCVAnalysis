# RCVAnalysis
Calculating features of rapid cyclic voltammetry traces

## Peak Metric Description
The following figure demonstrates a sample RCV scan with the peak metrics overlaid:
Performing RCV Analysis
In the ‘RCVAnalysis/build/bin’ subdirectory you can find a dynamic link library (RCVAnalysis.dll) which can be used to perform the analysis on raw RCV data. Run the peak finding function *findPeaksWrapper* to return indexes of the peaks followed by the feature calculation function *findFeaturesWrapper* to calculate the peak metrics demonstrated in section 1.1 (Peak Metric Description). These C++ wrapper functions take simple data types which make them more versatile with embedded systems.
NOTE: the raw RCV output needs to be formatted before passing it to RCVAnalysis.dll. xData is taken as the absolute values of the raw potential data and yData is taken as current in micro amps, converted using the following calibration equation:
  
*Current (μA) = 211.24 × Current(Potentiostat Output) + 0.455*
  
The HISENTS control application uses the RCVAnalysis.dll to analyse RCV data. The application outputs the metrics as a comma delimited file with the columns and rows detailed below. Note that the column and row headers are included here for reference only, these are not present in the HISENTS output files.



## RCV Analysis Methods
### Peak Finding Function
Function arguments for *findPeaksWrapper* are as follows:  
  
int **findPeaksWrapper**(double* xData,  
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

The HISENTS control application passes the following function argument values to *findPeaksWrapper*:
Function Argument | Passed Value
------------------|-------------
slopeThreshold | -1
ampThreshold | -Inf
smoothWidth | 5
peakGroup | 5
heightThreshold | 2.5
locationThreshold | 0.025

The method for finding peaks has been adapted from the ‘findpeaksxw.m’ script documented here: https://terpconnect.umd.edu/~toh/spectrum/PeakFindingandMeasurement.htm 
1. The RCV data is smoothed with a moving band of size *smoothWidth*
1. The derivative of the smoothed data is computed
1. The smoothed derivative is iterated in the positive direction, if a negative change is detected which has value greater than *slopeThreshold*, an array of size *peakGroup* is taken starting from the index at which the change was detected. The mean of the values from the peak group is stored as a peak candidate
1. Candidate peaks are eliminated if they have a height of less than *ampThreshold*
1. Candidate peaks are eliminated if they have an amplitude of value *heightThreshold* less than their positively indexed neighbour
1. Candidate peaks are eliminated if they have neighbours along the y-axis closer than value *locationThreshold*
1. The filtered peak indexes are returned for use with the feature calculation function

### Feature Calculation Function
Function arguments for *findFeaturesWrapper* are as follows:  
  
int **findFeaturesWrapper**(double* xData,  
double* yData,  
int* indexes,  
const unsigned int smoothWidth,  
double peakArea[2],  
double fullHeight[2],  
double FWHM[2],  
double xValsFullHeight[4],  
double xValsHalfHeight[4] )

Following smoothing and derivative calculation, steps 3-8 are repeated for each of the peak indexes identified using the peak finding function.
1. The RCV data is smoothed with a moving band of size *smoothWidth*
1. The derivative of the smoothed data is computed
1. Starting at the peak index, the smoothed derivative is iterated in the negative direction until a zero-crossing is identified. The y-value at this location is stored as the lower *fullHeight* amplitude (yValsFullHeight[0]).
1. Starting at the peak index, the smoothed derivative is iterated in the positive direction until a zero-crossing is identified. The y-value at this location is stored as the upper *fullHeight* amplitude (yValsFullHeight[1]).
1. *fullHeight* is calculated as: *((peakAmplitude - yValsFullHeight[0]) + (peakAmplitude - yValsFullHeight[1])) / 2*
1. *halfHeight* x-values are identified which correspond to y-locations either side of the peak with a y-value equal to *peakAmplitude - fullHeight / 2*
1. *FWHM* is calculated as the upper *halfHeight* x-value minus the lower *halfHeight* x-value
1. *peakArea* is calculated using the trapezoidal rule across the length of the peak (from the lower *fullHeight* x-value to the upper *fullHeight* x-value)
