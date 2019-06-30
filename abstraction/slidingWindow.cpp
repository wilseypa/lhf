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


// Find the minimum and maximum squared nearest neighbor distances within the window. These distances are computed
// periodically, for example, with the addition of every 50 new data points to the window. For every point in the
// window, compute and store the squared distance to its 2nd nearest neighbor (since every query point is within the
// window, the 1st nearest neighbor is the query point itself).
auto nnDistsWindow( vector<vector<double>> window ) {
    int k = 2;  // The number of near neighbors to search for.
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

        squaredNNdist.push_back(dists[1]);  // Store the squared 2nd nearest neighbor distance of the i-th point in the window.
    }

    // Clean up.
    delete [] nnIdx;
    delete [] dists;
    delete kdTree;
    annClose();

    // Find and return the minimum and maximum of all squared nearest neighbor distances within the window.
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

    int nnIndexInWindow = nnIdx[0];  // The index of the point in the window that is nearest to the incoming (query) point.
    double squaredNNDist = dists[0];  //  The squared distance to the point in the window that is nearest to the incoming (query) point.

    // Clean up.
    delete [] nnIdx;
    delete [] dists;
    delete kdTree;
    annClose();

    return tuple(nnIndexInWindow, squaredNNDist);  // Using the template argument deduction feature of C++17.
}

int main()
{
    unsigned int windowMaxSize{ 200 };  // The maximum number of points the sliding window may contain.
                            // It'll be a user input in future.

    unsigned int windowMinSize = windowMaxSize - 50;  // The minimum number of points that should be present in the window
                                             // to form meaningful topological features by PH computation.

    unsigned int numPointsAddedToWindow{ 0 };
    unsigned int nnDistCheckIntrvl{ 50 };

    // Variables to store the minimum and maximum of the squared nearest neighbor distances in the window.
    double minSqrdNNdist{ 0.0 };
    double maxSqrdNNdist{ 0.0 };

    int updateCounter = 0;


    vector<vector<double>> window;
    FILE *pFile;

    char buffer[1000];  // Define an arbitrary but large enough character array
                        // to store a single row of the data.

    pFile = fopen("abstraction/myFile.csv", "r");

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
                window.push_back(dataPoint);   // Add the new point to the back of the window.
                numPointsAddedToWindow++;
            }

            // Once window has more than the minimum number of points, apply a criterion for adding points in the window.
            else {
                if ( numPointsAddedToWindow % nnDistCheckIntrvl == 0 )  // If 50 new points have been added to the window:
                    tie(minSqrdNNdist, maxSqrdNNdist) = nnDistsWindow( window );  // compute the min and max of the squared
                                                                                  // nearest neighbor distances of the present
                                                                                  // set of points in the window.

                // Nearest neighbor search for an incoming point by ANN library. Using the structured binding feature of
                // C++17 to receive multiple values returned by the 'nnSearch' function.
                auto [nnIndex, sqrdDistToNearestRep] = nnSearch( dataPoint, window );


                // If the squared distance from the new incoming point to its nearest point in the window is
                // outside the range [minSqrdNNdist, maxSqrdNNdist]:
                if ( (sqrdDistToNearestRep < minSqrdNNdist) || (sqrdDistToNearestRep > maxSqrdNNdist) ) {
                    updateCounter++;
                    cout << updateCounter << ": " << minSqrdNNdist << ", " << maxSqrdNNdist << '\n';
                    window.push_back(dataPoint);  // add the incoming point to the back of the window.
                    numPointsAddedToWindow++;

                    // If the number of points in the window exceeds its maximum allowed size:
                    if ( window.size() > windowMaxSize ) {
                        vector<double> repToBeDeleted = window[0];
                        window.erase(window.begin());  // delete the point from the front of the window.
                    }
                }
                else {  // Discard the incoming point, and move its representative (nearest neighbor) to the back of the window.
                    auto representative = window.begin() + nnIndex;  // Iterator pointing to the representative of the incoming point.

                    // Move the representative to the back of the window. Since the majority of the incoming points is expected
                    // to find a representative than being added to the window, the 'move to back' operation is expected to be more
                    // frequent than the addition or deletion of points at either end of the window. Hence, the window is designed as
                    // a std::vector as opposed to std::deque. More information on the choice of std::vector can be found below.
                    // https://stackoverflow.com/questions/14579957/std-container-c-move-to-front
                    // https://stackoverflow.com/questions/20107756/push-existing-element-of-stddeque-to-the-front
                    rotate( representative, representative + 1, window.end() );  // https://stackoverflow.com/questions/23789498/moving-a-vector-element-to-the-back-of-the-vector
                }
            }

        }
        fclose (pFile);
    }
    return 0;
}
