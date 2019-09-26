#include <iostream>
#include <stdio.h>
#include <vector>
#include <string>
#include <iterator>
#include <sstream>
#include <algorithm>
#include <tuple>
#include <map>
#include <cmath>
#include <ANN/ANN.h>
#include "utils.hpp"
#include "basePipe.hpp"
#include "pipePacket.hpp"

// The next two functions are for splitting a string by a delimiter. They are
// used to split each row of the data by the appropriate delimiter.
// Source: https://stackoverflow.com/questions/236129/how-do-i-iterate-over-the-words-of-a-string
template<typename Out>
void split( const std::string &s, char delim, Out result ) {
    std::stringstream ss(s);
    std::string item;
    while ( std::getline(ss, item, delim) ) {
        *(result++) = item;
    }
}

std::vector<std::string> split( const std::string &s, char delim ) {
    std::vector<std::string> elems;
    split( s, delim, std::back_inserter(elems) );
    return elems;
}

// Convert a vector of strings to a vector of doubles. Each data point is
// represented as a vector of doubles.
// Source: https://stackoverflow.com/questions/20257582/convert-vectorstdstring-to-vectordouble
std::vector<double> stringVectorToDoubleVector( const std::vector<std::string>& stringVector ) {
    std::vector<double> doubleVector( stringVector.size() );
    std::transform( stringVector.begin(), stringVector.end(), doubleVector.begin(), [](const std::string& val)
                   {
                       return stod(val);
                   } );

    return doubleVector;
}

/*
// https://stackoverflow.com/questions/110157/how-to-retrieve-all-keys-or-values-from-a-stdmap-and-put-them-into-a-vector
std::vector<std::vector<double>> retrieveValuesFromMap( const std::map<int, std::vector<double>> &refWindow ) {
    std::vector<std::vector<double>> windowVals;
    windowVals.reserve( refWindow.size() );
    for (auto const &item : refWindow)
        windowVals.push_back(item.second);

    return windowVals;
}
*/

auto populateDistMatrix( pipePacket inData, std::vector<std::vector<double>> &refDistMat, std::vector<double> &refNNDists, double maxEpsilon ) {
    utils ut;
    for(unsigned int i = 0; i < inData.originalData.size(); i++) {
        std::vector<double> distsFromCurrVect;
        for(unsigned int j = 0; j < inData.originalData.size(); j++) {
            if (j < i) {
                distsFromCurrVect.push_back( refDistMat[j][i] );
            } else if (j > i) {
                auto dist = ut.vectors_distance( inData.originalData[i], inData.originalData[j] );
                if(dist < maxEpsilon)
                    inData.weights.insert( dist );
                refDistMat[i][j] = dist;
                distsFromCurrVect.push_back( dist );
            }
        }
        auto nnDistFromCurrVect = *std::min_element( distsFromCurrVect.begin(), distsFromCurrVect.end() );
        refNNDists.push_back( nnDistFromCurrVect );
    }

    inData.complex->setDistanceMatrix( refDistMat );
    inData.weights.insert( 0.0 );
    inData.weights.insert( maxEpsilon );

    auto nnDistsRange = std::minmax_element( refNNDists.begin(), refNNDists.end() );
    return std::tuple( *nnDistsRange.first, *nnDistsRange.second, inData );
}

void deleteFromDistMatrix( int indexToBeDeleted, std::vector<std::vector<double>> &refDistMat, std::vector<double> &refNNDists ) {
    refDistMat.erase( refDistMat.begin() + indexToBeDeleted );
    for(unsigned int i = 0; i < refDistMat.size(); i++)
        refDistMat[i].erase( refDistMat[i].begin() + indexToBeDeleted );

    refNNDists.erase( refNNDists.begin() + indexToBeDeleted );
    // auto nnDistsRange = std::minmax_element( refNNDists.begin(), refNNDists.end() );
    // return std::tuple( *nnDistsRange.first, *nnDistsRange.second );
}

auto addToDistMatrix( pipePacket inData, std::vector<std::vector<double>> &refDistMat, std::vector<double> &refNNDists, const std::vector<double> &refDistsFromNewPoint, double maxEpsilon ) {
    for(unsigned int i = 0; i < refDistMat.size(); i++) {
        refDistMat[i].push_back( refDistsFromNewPoint[i] );
        if( refDistsFromNewPoint[i] < maxEpsilon )
            inData.weights.insert( refDistsFromNewPoint[i] );
    }

    unsigned int newSize = refDistMat[0].size();
    std::vector<double> newRow(newSize, 0);
    refDistMat.push_back( newRow );

    inData.complex->setDistanceMatrix( refDistMat );

    auto nnDistFromNewPoint = *std::min_element( refDistsFromNewPoint.begin(), refDistsFromNewPoint.end() );
    refNNDists.push_back( nnDistFromNewPoint );

    auto nnDistsRange = std::minmax_element( refNNDists.begin(), refNNDists.end() );
    return std::tuple( *nnDistsRange.first, *nnDistsRange.second, inData );
}

/*
// Find the minimum and maximum squared nearest neighbor distances within the window. These distances are computed
// periodically, for example, with the addition of every 50 new data points to the window. For every point in the
// window, compute and store the squared distance to its 2nd nearest neighbor (since every query point is within the
// window, the 1st nearest neighbor is the query point itself).
auto nnDistsWindow( std::vector<std::vector<double>> window ) {
    int k = 2;  // The number of near neighbors to search for.
    int dim = window[0].size();  // Assuming all data points have the same dimension.
    double del = 0;  // At this time, we are using the exact nearest neighbor search.
    int nPts = window.size();  // The actual number of data points to search from.
    std::vector<double> squaredNNdist;  // A container to store the squared nearest neighbor distance
                                   // of each data point in the window.

    // Make the type of 'window' compatible with the type of an array of points as defined in the ANN library.
    // Source: https://stackoverflow.com/questions/4776164/c-vectorvectordouble-to-double
    std::vector<double*> ptrs;
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

    for(unsigned int i = 0; i < (unsigned)nPts; i++) {  // Iterate through each data point in the window.
        ANNpoint queryPt = &window[i][0];  // Each data point becomes a query point.
        kdTree->annkSearch(
                           queryPt,
                           k,
                           nnIdx,
                           dists,
                           del);

        squaredNNdist.push_back(dists[1]);  // Store the squared 2nd nearest neighbor distance of the i-th point in the window.
    }

    // Clean up.
    delete [] nnIdx;
    delete [] dists;
    delete kdTree;
    annClose();

    // Find and return the minimum and maximum of all squared nearest neighbor distances within the window.
    auto sqrdNNdistsRange = std::minmax_element( squaredNNdist.begin(), squaredNNdist.end() );
    return std::tuple( *sqrdNNdistsRange.first, *sqrdNNdistsRange.second );
}
*/

// Find the nearest neighbor of a new data point in the window. Return the index of the nearest neighbor
// and the nearest neighbor distance.
auto nnSearch( std::vector<double> &refNewPoint, std::vector<std::vector<double>> &refWindowValues ) {
    int k = refWindowValues.size();  // The number of near neighbors to search for.
    int dim = refNewPoint.size();
    double del = 0;  // At this time, we are using the exact nearest neighbor search.
    int nPts = refWindowValues.size();  // The actual number of data points to search from.

    // Make the type of 'window' compatible with the type of an array of points as defined in the ANN library.
    // Source: https://stackoverflow.com/questions/4776164/c-vectorvectordouble-to-double
    std::vector<double*> ptrs;
    for (auto &vec : refWindowValues)
        ptrs.push_back( vec.data() );

    ANNpointArray dataPts = ptrs.data();  // The data points to search from.

    ANNpoint queryPt = &refNewPoint[0];  // The query point which is the new incoming point from the stream.
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
                       del);

    // int nnIndexInWindow = nnIdx[0];  // The index of the point in the window that is nearest to the incoming (query) point.
    // double nnDistFromNewPoint = sqrt(dists[0]);  //  The distance to the point in the window that is nearest to the incoming (query) point.

    // Clean up.
    // delete [] nnIdx;
    // delete [] dists;
    delete kdTree;
    annClose();

    return std::tuple(nnIdx, dists);  // Using the template argument deduction feature of C++17.
}

int main()
{
    unsigned int windowMaxSize{ 200 };  // The maximum number of points the sliding window may contain.
                                        // It'll be a user input in future.

    unsigned int windowMinSize{ 150 };  // The minimum number of points that should be present in the window
                                        // to form meaningful topological features by PH computation.

    // unsigned int numPointsAddedToWindow{ 0 };
    // unsigned int nnDistCheckIntrvl{ 50 };

    // Variables to store the minimum and maximum of the squared nearest neighbor distances in the window.
    double minNNdist{ 0.0 };
    double maxNNdist{ 0.0 };

    std::vector<double> nnDists;
    nnDists.reserve( windowMaxSize );

    int addCounter{ 0 };
    int deleteCounter{ 0 };
    int lruCounter{ 0 };


    // std::map<int, std::vector<double>> window;
    std::vector<int> windowKeys;
    windowKeys.reserve( windowMaxSize );

    std::vector<std::vector<double>> windowValues;
    windowValues.reserve( windowMaxSize );

    std::vector<int> dynamicKeyContainer;
    dynamicKeyContainer.reserve( windowMaxSize );

    int key{ 0 };

    std::vector<std::vector<double>> distMatrix( windowMinSize, std::vector<double>(windowMinSize, 0) );
    distMatrix.reserve( windowMaxSize );

    std::map<std::string, std::string> args = { {"dimensions", "2"}, {"epsilon", "5"}, {"complexType", "indSimplexTree"} };
    // std::map<std::string, std::string> args = { {"dimensions", "2"}, {"epsilon", "5"}, {"complexType", "simplexTree"} };
    auto *wD = new pipePacket(args["complexType"], stod(args["epsilon"]), stoi(args["dimensions"]));

    double maxEpsilon = std::atof(args["epsilon"].c_str());

    std::cout << args["complexType"] << '\n';

    // auto *bp = new basePipe();
    // auto *cp = bp->newPipe( "neighGraph", args["complexType"] );


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
            std::vector<std::string> row = split( buffer, ',' );
            std::vector<double> dataPoint = stringVectorToDoubleVector(row);

            // std::copy( dataPoint.begin(), dataPoint.end(), std::ostream_iterator<double>(std::cout, " ") );
            // std::cout << '\n';

            // Initialize the window. Ensure that the sliding window always contains the minimum number of points.
            if ( windowKeys.size() < windowMinSize ) {
                // window.insert( std::pair<int, std::vector<double>>(label, dataPoint) );   // Add the new point to the back of the window.
                windowKeys.push_back(key);
                windowValues.push_back(dataPoint);
                dynamicKeyContainer.push_back(key);
                key++;
                // numPointsAddedToWindow++;
                if ( key == windowMinSize ) {
                    // std::vector<std::vector<double>> vectorsInWindow = retrieveValuesFromMap( window );
                    wD->originalData = windowValues;
                    pipePacket inData = *wD;
                    std::tie( minNNdist, maxNNdist, *wD ) = populateDistMatrix( *wD, distMatrix, nnDists, maxEpsilon );

                    auto *bp = new basePipe();
                    auto *cp = bp->newPipe( "neighGraph", args["complexType"] );

                    if(cp != 0 && cp->configPipe(args))
                        *wD = cp->runPipeWrapper(*wD);

                }
            }

            // Once window has more than the minimum number of points, apply a criterion for adding points in the window.
            else {
                // if ( numPointsAddedToWindow % nnDistCheckIntrvl == 0 )  // If 50 new points have been added to the window:
                    // std::tie(minSqrdNNdist, maxSqrdNNdist) = nnDistsWindow( window );  // compute the min and max of the squared
                                                                                  // nearest neighbor distances of the present
                                                                                  // set of points in the window.

                // Nearest neighbor search for an incoming point by ANN library. Using the structured binding feature of
                // C++17 to receive multiple values returned by the 'nnSearch' function.
                auto [repIndices, sqrdDistsToReps] = nnSearch( dataPoint, windowValues );
                double distToNearestRep = sqrt( sqrdDistsToReps[0] );

                // If the distance from the new incoming point to its nearest point in the window is
                // outside the range [minNNdist, maxNNdist]:
                if ( (distToNearestRep < minNNdist) || (distToNearestRep > maxNNdist) ) {
                    addCounter++;
                    std::cout << "AddCounter: " << addCounter << '\n';

                    // If the number of points in the window exceeds its maximum allowed size:
                    if ( windowKeys.size() == windowMaxSize ) {
                        deleteCounter++;
                        std::cout << "DeleteCounter: " << deleteCounter << '\n';
                        int keyToBeDeleted = dynamicKeyContainer[0];
                        dynamicKeyContainer.erase( dynamicKeyContainer.begin() );  // delete the point from the front of the window.

                        // https://thispointer.com/c-how-to-find-an-element-in-vector-and-get-its-index/
                        std::vector<int>::iterator iter = std::find( windowKeys.begin(), windowKeys.end(), keyToBeDeleted );
                        int indexToBeDeleted = std::distance( windowKeys.begin(), iter );

                        windowKeys.erase( windowKeys.begin() + indexToBeDeleted );
                        windowValues.erase( windowValues.begin() + indexToBeDeleted );
                        deleteFromDistMatrix( indexToBeDeleted, distMatrix, nnDists );

                        //std::set<unsigned> deletionSet = {indexToBeDeleted};
                        //wd->complex->deletion(deletionSet);


                    }

                    windowKeys.push_back(key);  // add the incoming point to the back of the window.
                    windowValues.push_back(dataPoint);
                    dynamicKeyContainer.push_back(key);
                    key++;
                    // numPointsAddedToWindow++;

                    std::vector<double> distsFromNewPoint;
                    unsigned int numReps = windowKeys.size() - 1;
                    for(unsigned int i = 0; i < numReps; i++) {
                        // https://stackoverflow.com/questions/3909784/how-do-i-find-a-particular-value-in-an-array-and-return-its-index
                        int windowIndex = std::distance(repIndices, std::find(repIndices, repIndices + numReps, i));
                        distsFromNewPoint.push_back( sqrt( sqrdDistsToReps[windowIndex] ) );
                    }

                    wD->originalData.push_back( dataPoint );
                    std::tie(minNNdist, maxNNdist, *wD) = addToDistMatrix( *wD, distMatrix, nnDists, distsFromNewPoint, maxEpsilon );
                    wD->complex->insert( dataPoint );

                }
                else {  // Discard the incoming point, and move its representative (nearest neighbor) to the back of the window.
                    int nnRepIndex = repIndices[0];
                    int nnRepKey = windowKeys[(unsigned)nnRepIndex];

                    std::vector<int>::iterator itRep = std::find( dynamicKeyContainer.begin(), dynamicKeyContainer.end(), nnRepKey ); // Iterator pointing to the representative of the incoming point.

                    // Move the representative to the back of the window. Since the majority of the incoming points is expected
                    // to find a representative than being added to the window, the 'move to back' operation is expected to be more
                    // frequent than the addition or deletion of points at either end of the window. Hence, the window is designed as
                    // a std::vector as opposed to std::deque. More information on the choice of std::vector can be found below.
                    // https://stackoverflow.com/questions/14579957/std-container-c-move-to-front
                    // https://stackoverflow.com/questions/20107756/push-existing-element-of-stddeque-to-the-front
                    std::rotate( itRep, itRep + 1, dynamicKeyContainer.end() );  // https://stackoverflow.com/questions/23789498/moving-a-vector-element-to-the-back-of-the-vector
                    lruCounter++;
                    // std::cout << "LRU Counter: " << lruCounter << '\n';
                }
            }

        }
        fclose (pFile);
    }

    return 0;
}
