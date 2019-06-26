#include <iostream>
#include <stdio.h>
#include <vector>
#include <string>
#include <iterator>
#include <sstream>
#include <algorithm>
#include <tuple>
#include <ANN/ANN.h>

using namespace std;

// The next two functions are for splitting a string by a delimiter. They are
// used to split each row of the data by the appropriate delimiter.
// Source: https://stackoverflow.com/questions/236129/how-do-i-iterate-over-the-words-of-a-string
template<typename Out>
void split( const string &s, char delim, Out result ) {
    stringstream ss(s);
    string item;
    while ( getline(ss, item, delim) ) {
        *(result++) = item;
    }
}

vector<string> split( const string &s, char delim ) {
    vector<string> elems;
    split( s, delim, back_inserter(elems) );
    return elems;
}

// Convert a vector of strings to a vector of doubles. Each data point is
// represented as a vector of doubles.
// Source: https://stackoverflow.com/questions/20257582/convert-vectorstdstring-to-vectordouble
vector<double> stringVectorToDoubleVector( const vector<string>& stringVector ) {
    vector<double> doubleVector( stringVector.size() );
    transform( stringVector.begin(), stringVector.end(), doubleVector.begin(), [](const string& val)
                   {
                       return stod(val);
                   } );

    return doubleVector;
}


auto nnDistsWindow( vector<vector<double>> window ) {
    int k = 1;  // The number of near neighbors to search for.
    int dim = window[0].size();  // Assuming all data points have the same dimension.
    double eps = 0;  // At this time, we are using the exact nearest neighbor search.
    int nPts = window.size();  // The actual number of data points to search from.
    vector<double> squaredNNdist;  // A container to store the squared nearest neighbor distance
                                   // of each data point in the window.

    // Make the type of 'window' compatible with the type of an array of points as defined in the ANN library.
    // Source: https://stackoverflow.com/questions/4776164/c-vectorvectordouble-to-double
    vector<double*> ptrs;
    for (auto& vec : window)
        ptrs.push_back(vec.data());

    ANNpointArray dataPts = ptrs.data();  // The data points to search from.

    ANNidxArray nnIdx = new ANNidx[k];  // Near neighbor indices (to be returned by the search).
    ANNdistArray dists = new ANNdist[k];  // Squared near neighbor distances (to be returned by the search).

    // Build the search structure.
    ANNkd_tree* kdTree = new ANNkd_tree(
                            dataPts,
                            nPts,
                            dim);

    for (unsigned int i = 0; i < (unsigned)nPts; i++) {  // Iterate through each data point in the window.
        ANNpoint queryPt = &window[i][0];  // Each data point becomes a query point.
        kdTree->annkSearch(
                           queryPt,
                           k,
                           nnIdx,
                           dists,
                           eps);

        squaredNNdist.push_back(dists[0]);  // Store the squared nearest neighbor distance of the i-th point in the window.
    }

    auto sqrdNNdistsRange = minmax_element( squaredNNdist.begin(), squaredNNdist.end() );
    return tuple( *sqrdNNdistsRange.first, *sqrdNNdistsRange.second );
}

// Find the nearest neighbor of a new data point in the window. Return the index of the nearest neighbor
// and the squared nearest neighbor distance.
auto nnSearch( vector<double> newPoint, vector<vector<double>> window ) {
    int k = 1;  // The number of near neighbors to search for.
    int dim = newPoint.size();
    double eps = 0;  // At this time, we are using the exact nearest neighbor search.
    int nPts = window.size();  // The actual number of data points to search from.

    // Make the type of 'window' compatible with the type of an array of points as defined in the ANN library.
    // Source: https://stackoverflow.com/questions/4776164/c-vectorvectordouble-to-double
    vector<double*> ptrs;
    for (auto& vec : window)
        ptrs.push_back(vec.data());

    ANNpointArray dataPts = ptrs.data();  // The data points to search from.

    ANNpoint queryPt = &newPoint[0];  // The query point which is the new incoming point from the stream.
    ANNidxArray nnIdx = new ANNidx[k];  // Near neighbor indices (to be returned by the search).
    ANNdistArray dists = new ANNdist[k];  // Squared near neighbor distances (to be returned by the search).

    // Build the search structure.
    ANNkd_tree* kdTree = new ANNkd_tree(
                            dataPts,
                            nPts,
                            dim);

    kdTree->annkSearch(
                       queryPt,
                       k,
                       nnIdx,
                       dists,
                       eps);

    int nnIndexinWindow = nnIdx[0];
    double squaredNNDist = dists[0];

    // Clean up.
    delete [] nnIdx;
    delete [] dists;
    delete kdTree;
    annClose();

    return tuple(nnIndexinWindow, squaredNNDist);  // Using the template argument deduction feature of C++17.
}

int main()
{
    unsigned int windowMaxSize{ 200 };  // The maximum number of points the sliding window may contain.
                            // It'll be a user input in future.

    unsigned int windowMinSize = windowMaxSize - 50;  // The minimum number of points that should be present in the window
                                             // to form meaningful topological features by PH computation.

    unsigned int numPointsAddedToWindow{ 0 };
    unsigned int nnDistCheckIntrvl{ 50 };

    vector<vector<double>> window;
    FILE *pFile;

    char buffer[1000];  // Define an arbitrary but large enough character array
                        // to store a single row of the data.

    pFile = fopen("myFile.csv", "r");

    if (pFile == NULL)
        perror ("Error opening file");
    else {
        while( !feof(pFile) ) {
            if ( fgets(buffer, 1000, pFile) == NULL )  // Read one row at a time and store in 'buffer'.
                break;

            // Split the row by the appropriate delimiter to store it as a vector of doubles.
            vector<string> row = split( buffer, ',' );
            vector<double> dataPoint = stringVectorToDoubleVector(row);

            // copy( dataPoint.begin(), dataPoint.end(), ostream_iterator<double>(cout, " ") );
            // cout << '\n';

            // Initialize the window. Ensure that the sliding window always contains the minimum number of points.
            if ( window.size() < windowMinSize ) {
                window.push_back(dataPoint);
                numPointsAddedToWindow++;
            }

            // Once window has more than the minimum number of points, apply a criterion for adding points in the window.
            else {
                if ( numPointsAddedToWindow % nnDistCheckIntrvl == 0 )
                    auto [minSqrdNNdist, maxSqrdNNdist] = nnDistsWindow( window );

                // Nearest neighbor search for an incoming point by ANN library. Using the structured binding feature of
                // C++17 to receive multiple values returned by the 'nnSearch' function.
                auto [nnIndex, sqrdDistToNearestRep] = nnSearch( dataPoint, window );
                if ( (sqrdDistToNearestRep < minSqrdNNdist) || (sqrdDistToNearestRep > maxSqrdNNdist) ) {
                    window.push_back(dataPoint);
                    numPointsAddedToWindow++;
                    if ( window.size() > windowMaxSize )
                }
            }


        }
        fclose (pFile);
    }
    return 0;
}
