/*
 * streamVR hpp + cpp extend the basePipe class for
 *
 */

#include <string>
#include <cstring>   //added for memory copy
#include <iostream>
#include <fstream>
#include <iterator>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <numeric>
#include <functional>
#include "slidingWindow.hpp"
#include "readInput.hpp"


int winp = 0;
int pointCounter = 1;
std::vector<std::vector<double>> distMatrix;

std::vector<int> avgfirst;
std::vector<double> avgsecond;

std::vector<int> numPointfirst;
std::vector<int> numPointsecond;

std::vector<int> maxKeyfirst;
std::vector<int> maxKeysecond;


//slidingWindow::EvalParams* slidingWindow::defaultVals = new EvalParams{ 10, 0, 0 };
slidingWindow::EvalParams* slidingWindow::defaultVals = new EvalParams{ WINDOWSIZE, 0, 0 };

slidingWindow::EvalParams* slidingWindow::totalDefaultVals=new EvalParams{WINDOWSIZE*WINDOWNUMS+1,0,0};

//initialize for the sumDefaultVals
//slidingWindow::EvalParams* slidingWindow::sumDefaultVals0 = new EvalParams{ 10, 0, 0 };  //initalization and reference
slidingWindow::EvalParams* slidingWindow::sumDefaultVals0 = new EvalParams{ WINDOWSIZE, 0, 0 }; 

slidingWindow::EvalParams* slidingWindow::sumDefaultValsT[WINDOWNUMS+1];  //reference

pipePacket* slidingWindow::pPack;
pipePacket* slidingWindow::sumpPackT[WINDOWNUMS+1];
pipePacket* slidingWindow::totalpPack;

//std::vector<std::vector<double>> slidingWindow::sumoriginalData[WINDOWNUMS+1];

std::vector<std::vector<std::vector<double>>> slidingWindow::sumoriginalData;

// basePipe constructor
slidingWindow::slidingWindow()
{
	slidingWindow::sumoriginalData.resize(WINDOWNUMS+1);
	for (int k=0; k<WINDOWNUMS+1; k++)
	{
	
	//slidingWindow::sumDefaultVals0 = new EvalParams{ 10, 0, 0 };
	slidingWindow::sumDefaultVals0 = new EvalParams{ WINDOWSIZE, 0, 0 };
	
    slidingWindow::sumDefaultValsT[k]=slidingWindow::sumDefaultVals0;
    
    
    }
	
    pipeType = "SlidingWindow";
    return;
}

void slidingWindow::deleteNNstats()
{
    defaultVals->keyToBeDeleted = defaultVals->windowKeys[defaultVals->indexToBeDeleted];
    
    std::cout << "keyToBeDeleted: " << defaultVals->keyToBeDeleted <<" ";

    defaultVals->windowKeys.erase( defaultVals->windowKeys.begin() + defaultVals->indexToBeDeleted );
  
    defaultVals->partitionLabels.erase( defaultVals->partitionLabels.begin() + defaultVals->indexToBeDeleted );
    
    

    defaultVals->nnIndices.erase( defaultVals->nnIndices.begin() + defaultVals->indexToBeDeleted );

    defaultVals->nnDistToBeDeleted = defaultVals->nnDists[defaultVals->indexToBeDeleted];
    defaultVals->nnDists.erase( defaultVals->nnDists.begin() + defaultVals->indexToBeDeleted );

    defaultVals->distsFromCurrVec.erase( defaultVals->distsFromCurrVec.begin() + defaultVals->indexToBeDeleted );

    defaultVals->numPointsPartn[defaultVals->labelToBeDeleted] = defaultVals->numPointsPartn[defaultVals->labelToBeDeleted] - 1;
    
   
    // If there are no more points left in the partition from which the deletion took place:
    if (defaultVals->numPointsPartn[defaultVals->labelToBeDeleted] == 0)
    {
        defaultVals->avgNNDistPartitions.erase(defaultVals->labelToBeDeleted);
        defaultVals->numPointsPartn.erase(defaultVals->labelToBeDeleted);
        defaultVals->maxKeys.erase(defaultVals->labelToBeDeleted);
    }

    return;

}

void slidingWindow::updateStats()
{
    // Delete the corresponding row from the distance matrix.
    // pPack->complex->distMatrix->erase(pPack->complex->distMatrix->begin() + defaultVals->indexToBeDeleted);
    distMatrix.erase(distMatrix.begin() + defaultVals->indexToBeDeleted);


    double sumOldNNdists = 0.0;
    double sumNewNNdists = 0.0;

    double sumDelOldNNdists = 0.0;
    double sumDelNewdists = 0.0;

    double sumAddOldNNdists = 0.0;
    double sumAddNewdists = 0.0;

    int twoPointPartition = 0;

    std::vector<double> distsFromCurrVecTP;  // A vector to store the distances from the current vector to the ones in the target partition.
    std::vector<int> tpIndices; // A vector to store the positions of the existing members of the target partition.

    // Add a new column and row to the end of the upper triangular distance matrix.
    for(unsigned int i = 0; i < distMatrix.size(); i++)
    {
        // Delete the corresponding entry from each row.
        // pPack->complex->distMatrix->at(i).erase( pPack->complex->distMatrix->at(i).begin() + defaultVals->indexToBeDeleted );
        distMatrix[i].erase( distMatrix[i].begin() + defaultVals->indexToBeDeleted );

        // Add the new distance value to the end of each row.
        // pPack->complex->distMatrix->at(i).push_back( defaultVals->distsFromCurrVec[i] );
        distMatrix[i].push_back( defaultVals->distsFromCurrVec[i] );

        // Update NN statistics for only those partitions from which the point was deleted or to which the new point is to be added.
        // Case 1: The i-th point belongs to the partition the last point was deleted from, but not to the partition the new point is
        // to be added to. And, the point that was deleted was the nearest neighbor of the i-th point.
        if ( defaultVals->partitionLabels[i] == defaultVals->labelToBeDeleted
                &&  defaultVals->partitionLabels[i] != defaultVals->targetPartition
                && defaultVals->nnIndices[i] == defaultVals->indexToBeDeleted )
        {

            // If only one point remains in the partition from which the last point was deleted:
            if ( defaultVals->numPointsPartn[defaultVals->labelToBeDeleted] == 1 )   // labelToBeDeketed is the label of the partition that has the last point be deleted
            {
                defaultVals->avgNNDistPartitions[defaultVals->labelToBeDeleted] = -1;
                defaultVals->nnIndices[i] = -1;
                defaultVals->nnDists[i] = -1;

            }
            else
            {
                sumDelOldNNdists = sumDelOldNNdists + defaultVals->nnDists[i];

                std::vector<double> memberDistsFromVect;
                std::vector<int> memberIndices;

                for(unsigned int j = 0; j < distMatrix.size(); j++)
                {
                    if ( defaultVals->partitionLabels[j] == defaultVals->labelToBeDeleted )
                    {
                        if (j < i)
                        {
                            memberDistsFromVect.push_back( distMatrix[j][i] );
                            memberIndices.push_back(j);
                        }
                        else if (j > i)
                        {
                            memberDistsFromVect.push_back( distMatrix[i][j] );
                            memberIndices.push_back(j);
                        }
                    }
                }
                auto tempIdx = std::min_element(memberDistsFromVect.begin(), memberDistsFromVect.end()) - memberDistsFromVect.begin();
                defaultVals->nnIndices[i] = memberIndices[tempIdx];

                auto newNNdistFromVect = *std::min_element( memberDistsFromVect.begin(), memberDistsFromVect.end() );
                defaultVals->nnDists[i] = newNNdistFromVect;
                sumDelNewdists = sumDelNewdists + newNNdistFromVect;
            }

        }

        // Case 2: The i-th point belongs to the partition the new point is to be added to, but not to the partition the last point
        // was deleted from. In this case, is it possible that the point that was deleted was the nearest neighbor of the i-th point?
        else if ( defaultVals->partitionLabels[i] != defaultVals->labelToBeDeleted && defaultVals->partitionLabels[i] == defaultVals->targetPartition )
        {

            // If, previously, there was only one point in the target partition:
            if ( defaultVals->numPointsPartn[defaultVals->targetPartition] == 1 )
            {

                // Update the NN statistics for the "new" two-point partition.
                defaultVals->nnIndices[i] = defaultVals->windowMaxSize - 1;
                defaultVals->nnIndices.push_back(i);

                defaultVals->nnDists[i] = defaultVals->distsFromCurrVec[i];
                defaultVals->nnDists.push_back(defaultVals->nnDists[i]);

                defaultVals->avgNNDistPartitions[defaultVals->targetPartition] = defaultVals->nnDists[i];
                defaultVals->numPointsPartn[defaultVals->targetPartition] = 2;
                twoPointPartition = 1;
            }

            else
            {

                distsFromCurrVecTP.push_back( defaultVals->distsFromCurrVec[i] );
                tpIndices.push_back(i);

                if ( defaultVals->distsFromCurrVec[i] < defaultVals->nnDists[i] )
                {
                    sumAddOldNNdists = sumAddOldNNdists + defaultVals->nnDists[i];

                    defaultVals->nnDists[i] = defaultVals->distsFromCurrVec[i];
                    defaultVals->nnIndices[i] = defaultVals->windowMaxSize - 1;

                    sumAddNewdists = sumAddNewdists + defaultVals->nnDists[i];
                }
            }

        }

        // Case 3: The partition the last point was deleted from is the same as the target partition of the new point to be added to,
        // and the i-th point belongs to that partition (i.e. at least one point remained in the partition after the deletion and before
        // the addition).
        else if (defaultVals->labelToBeDeleted == defaultVals->targetPartition && defaultVals->partitionLabels[i] == defaultVals->targetPartition)
        {

            // If only one point remained in the partition after the deletion and before the addition:
            if ( defaultVals->numPointsPartn[defaultVals->targetPartition] == 1 )
            {

                // Update the NN statistics for the "new" two-point partition.
                defaultVals->nnIndices[i] = defaultVals->windowMaxSize - 1;
                defaultVals->nnIndices.push_back(i);

                defaultVals->nnDists[i] = defaultVals->distsFromCurrVec[i];
                defaultVals->nnDists.push_back(defaultVals->nnDists[i]);

                defaultVals->avgNNDistPartitions[defaultVals->targetPartition] = defaultVals->nnDists[i];
                defaultVals->numPointsPartn[defaultVals->targetPartition] = 2;
                twoPointPartition = 1;
            }

            else if ( defaultVals->distsFromCurrVec[i] < defaultVals->nnDists[i] )
            {
                distsFromCurrVecTP.push_back( defaultVals->distsFromCurrVec[i] );
                tpIndices.push_back(i);

                sumOldNNdists = sumOldNNdists + defaultVals->nnDists[i];

                defaultVals->nnDists[i] = defaultVals->distsFromCurrVec[i];
                defaultVals->nnIndices[i] = defaultVals->windowMaxSize - 1;

                sumNewNNdists = sumNewNNdists + defaultVals->nnDists[i];
            }

            else if ( defaultVals->nnIndices[i] == defaultVals->indexToBeDeleted )
            {
                distsFromCurrVecTP.push_back( defaultVals->distsFromCurrVec[i] );
                tpIndices.push_back(i);

                sumOldNNdists = sumOldNNdists + defaultVals->nnDists[i];

                std::vector<double> memberDistsFromVect;
                std::vector<int> memberIndices;

                for(unsigned int j = 0; j < defaultVals->windowMaxSize; j++)
                {
                    if ( defaultVals->partitionLabels[j] == defaultVals->labelToBeDeleted )
                    {
                        if (j < i)
                        {
                            memberDistsFromVect.push_back( distMatrix[j][i] );
                            memberIndices.push_back(j);
                        }
                        else if (j > i)
                        {
                            memberDistsFromVect.push_back( distMatrix[i][j] );
                            memberIndices.push_back(j);
                        }
                    }
                }
                auto tempIdx = std::min_element(memberDistsFromVect.begin(), memberDistsFromVect.end()) - memberDistsFromVect.begin();
                defaultVals->nnIndices[i] = memberIndices[tempIdx];

                auto newNNdistFromVect = *std::min_element( memberDistsFromVect.begin(), memberDistsFromVect.end() );
                defaultVals->nnDists[i] = newNNdistFromVect;
                sumNewNNdists = sumNewNNdists + newNNdistFromVect;
            }

            else
            {
                distsFromCurrVecTP.push_back( defaultVals->distsFromCurrVec[i] );
                tpIndices.push_back(i);
            }
        }

    }

    std::vector<double> distMatLastRow(defaultVals->windowMaxSize);  // The last row of the upper triangular distance matrix is a vector of 0s.
    // std::vector<double> distMatLastRow = defaultVals->distsFromCurrVec;
    // distMatLastRow.push_back(0);
    //pPack->complex->distMatrix->push_back( distMatLastRow );
    distMatrix.push_back( distMatLastRow );

    pPack->complex->setDistanceMatrix(&distMatrix);


    // Update the average NN distance of the partition from which the last point was deleted and of the one to which the new point
    // is being added.
    // If the new point is being added to a partition that does not exist in the window (after the last deletion):
    if ( defaultVals->avgNNDistPartitions.count( defaultVals->targetPartition ) == 0 )
    {
        defaultVals->avgNNDistPartitions[defaultVals->targetPartition] = -1;
        defaultVals->nnIndices.push_back(-1);
        defaultVals->nnDists.push_back(-1);
        defaultVals->numPointsPartn[defaultVals->targetPartition] = 1;

    }
    if ( defaultVals->labelToBeDeleted != defaultVals->targetPartition )
    {

        // Update the average NN distance of the partition from which the last point was deleted.
        if ( defaultVals->numPointsPartn[defaultVals->labelToBeDeleted] > 1 )
        {

            int n = defaultVals->numPointsPartn[defaultVals->labelToBeDeleted];
            double avgDel = defaultVals->avgNNDistPartitions[defaultVals->labelToBeDeleted];
            defaultVals->avgNNDistPartitions[defaultVals->labelToBeDeleted] = ( (n+1)*avgDel - defaultVals->nnDistToBeDeleted - sumDelOldNNdists + sumDelNewdists ) / n;
        }

        // Update the average NN distance of the partition to which the new point is being added.
        if ( defaultVals->numPointsPartn[defaultVals->targetPartition] >= 2 && twoPointPartition == 0 )
        {
            auto tempNNidx = std::min_element(distsFromCurrVecTP.begin(), distsFromCurrVecTP.end()) - distsFromCurrVecTP.begin();
            defaultVals->nnIndices.push_back( tpIndices[tempNNidx] );

            auto nnDistFromNewVect = *std::min_element( distsFromCurrVecTP.begin(), distsFromCurrVecTP.end() );
            defaultVals->nnDists.push_back( nnDistFromNewVect );

            int n = defaultVals->numPointsPartn[defaultVals->targetPartition];
            double avgAdd = defaultVals->avgNNDistPartitions[defaultVals->targetPartition];
            defaultVals->avgNNDistPartitions[defaultVals->targetPartition] = ( n*avgAdd + nnDistFromNewVect - sumAddOldNNdists + sumAddNewdists ) / (n+1);

            defaultVals->numPointsPartn[defaultVals->targetPartition] = n + 1;

        }
    }
    else if ( defaultVals->labelToBeDeleted == defaultVals->targetPartition
              && defaultVals->numPointsPartn[defaultVals->targetPartition] >= 2
              && twoPointPartition == 0 )
    {
        auto tempNNidx = std::min_element(distsFromCurrVecTP.begin(), distsFromCurrVecTP.end()) - distsFromCurrVecTP.begin();
        defaultVals->nnIndices.push_back( tpIndices[tempNNidx] );

        auto nnDistFromNewVect = *std::min_element( distsFromCurrVecTP.begin(), distsFromCurrVecTP.end() );
        defaultVals->nnDists.push_back( nnDistFromNewVect );

        int n = defaultVals->numPointsPartn[defaultVals->targetPartition] + 1;
        double avgNNd = defaultVals->avgNNDistPartitions[defaultVals->targetPartition];

        defaultVals->avgNNDistPartitions[defaultVals->targetPartition] = ( n*avgNNd - defaultVals->nnDistToBeDeleted - sumOldNNdists + sumNewNNdists + nnDistFromNewVect ) / n;
        defaultVals->numPointsPartn[defaultVals->targetPartition] = n;
    }

    return;
}

bool slidingWindow::nnBasedEvaluator(std::vector<double>& currentVector, std::vector<std::vector<double>>& windowValues)
{
    utils ut;                                            // utils.cpp line 15
    float f1{ 4 };
    float f2{ 1 };
    float f3{ 5 };
    float f4{ 3 };
    defaultVals->distsFromCurrVec.clear();            //sildingWindow.hpp line 32: std::vector<double> distsFromCurrVec

    // Compute the distances from the current vector to the existing ones in the window.
    defaultVals->distsFromCurrVec = ut.nearestNeighbors(currentVector, windowValues);                     //utils.cpp line 558
		  
																   //std::vector<double> utils::nearestNeighbors(std::vector<double>& point, std::vector<std::vector<double>>& pointcloud)
																    //{ ////based on random projection, x is current point being examined, n is number of centroids/facilities
																       //utils ut;
																       //std::vector<double> retVal;
																	   ////Get sq distances for each point
																       //for(auto currentPoint : pointcloud)
																        //{
																	     //retVal.push_back(ut.vectors_distance(point, currentPoint));
																        //}
																       //return retVal;
																     //}                                                                                                 
                                                                                                          
    // Find the distance from the current vector to its nearest neighbor in the window.
    auto nnDistCurrVec = *std::min_element( defaultVals->distsFromCurrVec.begin(), defaultVals->distsFromCurrVec.end() );
    
    //csyoutDebug("distsFromCurrVec 384", defaultVals->distsFromCurrVec,PRCDEBUGGER);


    std::cout << " 384 nnDistCurrVec = " << nnDistCurrVec << std::endl;

    if (defaultVals->avgNNDistPartitions.size() == 1)  // If the window is 'pure':
    {

        if (nnDistCurrVec == 0) {

            std::cout << "pointCounter = " << pointCounter << '\n';
            for (auto const& pair: defaultVals->avgNNDistPartitions) {
                std::cout << "{" << pair.first << ": " << pair.second << "}\n";
            }
            std::cout << "1 Skipped ==============================================================================================" << '\n';

            return false;
        }

        // Find the average nearest neighbor distance in the single 'partition' in the window.
        int currentLabel = defaultVals->partitionLabels[0];
        // std::cout << "currentLabel: " << currentLabel << '\n';

        auto avgNNDistSinglePartition = defaultVals->avgNNDistPartitions[currentLabel];

        std::cout << "avgNNDistSinglePartition = " << avgNNDistSinglePartition << std::endl;

        if (avgNNDistSinglePartition <= f2 && nnDistCurrVec <= 1.5) {

            std::cout << "pointCounter = " << pointCounter << '\n';
            for (auto const& pair: defaultVals->avgNNDistPartitions) {
                std::cout << "{" << pair.first << ": " << pair.second << "}\n";
            }
            std::cout << "2 Skipped ==============================================================================================" << '\n';

            return false;
        }


        if (avgNNDistSinglePartition == 0 || nnDistCurrVec / avgNNDistSinglePartition > f1)
        {
            // In this case, the oldest point will be deleted from the window.
            defaultVals->labelToBeDeleted = defaultVals->partitionLabels[0];
            defaultVals->indexToBeDeleted = 0;

            // Since the window was 'pure', the partition label of the new point to be added is (existing label + 1).
            defaultVals->targetPartition = defaultVals->labelToBeDeleted + 1;

            deleteNNstats();

            // Delete the vector from the front of the sliding window.
            windowValues.erase( windowValues.begin() );

            // Update the distance matrix and the NN statistics.
            updateStats();

            std::cout << "1. Point Added" << '\n';

            std::cout << "pointCounter = " << pointCounter << '\n';
            for (auto const& pair: defaultVals->avgNNDistPartitions) {
                std::cout << "{" << pair.first << ": " << pair.second << "}\n";
            }
            std::cout << "==============================================================================================" << '\n';

            return true;
        }

        std::cout << "pointCounter = " << pointCounter << '\n';
        for (auto const& pair: defaultVals->avgNNDistPartitions)
        {
            std::cout << "{" << pair.first << ": " << pair.second << "}  ";
        }


        std::cout << "\nSkipped ==============================================================================================" << '\n';

    }

    else {   // If the window is NOT "pure":


        for (auto const& pair: defaultVals->avgNNDistPartitions)
        {
            std::cout << "{" << pair.first << ": " << pair.second << "}  ";
        }

        // Determine the partition membership of the current vector. In particular, check if the current vector can be assigned to
        // its nearest partition. If it cannot be assigned to its nearest partition, create a new partition with only the current vector.

        // Find the partition that is nearest to the current vector.
        auto nearestPartnIdx = std::min_element(defaultVals->distsFromCurrVec.begin(), defaultVals->distsFromCurrVec.end()) - defaultVals->distsFromCurrVec.begin();
        int nearestPartition = defaultVals->partitionLabels[nearestPartnIdx];

        // Find the max. of the existing partition labels.
        int maxLabel = *std::max_element( defaultVals->partitionLabels.begin(), defaultVals->partitionLabels.end() );

        // If the (nearest neighbor) distance from the current vector to any partition is 0, assign the current vector to that partition.
        if ( nnDistCurrVec == 0 ) {
           defaultVals->targetPartition = nearestPartition;
        }

        // If the minimum distance from the current vector to existing partitions is higher than f3, assign a new partition to the current vector.
        else if ( nnDistCurrVec > f3 ) {
            defaultVals->targetPartition = maxLabel + 1;
        }
        else if ( defaultVals->avgNNDistPartitions[nearestPartition] == -1 && nnDistCurrVec <= f4 ) {
            defaultVals->targetPartition = nearestPartition;
        }
        else if ( defaultVals->avgNNDistPartitions[nearestPartition] <= f2 && nnDistCurrVec <= 1.5 ) {
            defaultVals->targetPartition = nearestPartition;
        }
        else if ( defaultVals->avgNNDistPartitions[nearestPartition] > 0 && nnDistCurrVec / defaultVals->avgNNDistPartitions[nearestPartition] <= f1 ) {
            defaultVals->targetPartition = nearestPartition;
        }
        else {
            defaultVals->targetPartition = maxLabel + 1;
        }

        std::cout << "targetPartition = " << defaultVals->targetPartition << '\n';

        // Determine which point would be deleted from the sliding window.
        // A partition is considered outdated if it did not receive any new point for more than the last 20 insertions.
        int timeToBeOutdated{ 10 };

        // First check if there are outdated partitions. If there are, find the smallest one of them.
        int smallestOutdated{ -1 };   // This will store the label of the smallest outdated partition.
        int numPtsSmallestOutdated{ -1 };  // This will store the number of points in the smallest outdated partition.

        for(auto op : defaultVals->maxKeys) {
            if ( (defaultVals->key - op.second) > timeToBeOutdated ) {   // If op is outdated:
                if ( numPtsSmallestOutdated == -1 ) {   // If this is the first time an outdated partition is encountered:
                    numPtsSmallestOutdated = defaultVals->numPointsPartn[op.first];
                    smallestOutdated = op.first;
                }
                else if ( defaultVals->numPointsPartn[op.first] < numPtsSmallestOutdated ) {   // If a smaller outdated partition is encountered:
                    // Update our records for the smallest outdated partition.
                    numPtsSmallestOutdated = defaultVals->numPointsPartn[op.first];
                    smallestOutdated = op.first;
                }
            }
        }

        if (smallestOutdated != -1) {   // If there is (are) outdated partition(s):

            // Delete the oldest point of the smallest outdated partition (and its associated statistics) from the sliding window.
            defaultVals->labelToBeDeleted = smallestOutdated;
            defaultVals->indexToBeDeleted = std::find( defaultVals->partitionLabels.begin(), defaultVals->partitionLabels.end(), smallestOutdated ) - defaultVals->partitionLabels.begin();

            deleteNNstats();

            windowValues.erase( windowValues.begin() + defaultVals->indexToBeDeleted );

            // Update the distance matrix and the NN statistics.
            updateStats();

            std::cout << "2. Point Added" << '\n';

            std::cout << "pointCounter = " << pointCounter << '\n';
            for (auto const& pair: defaultVals->avgNNDistPartitions) {
                std::cout << "{" << pair.first << ": " << pair.second << "}\n";
            }
            std::cout << "==============================================================================================" << '\n';

        }

        else {   // There is no outdated partition in the window:
            if ( defaultVals->targetPartition != nearestPartition ) {   // If the current vector was assigned a new partition:

                // In this case, the oldest point will be deleted from the window.
                defaultVals->labelToBeDeleted = defaultVals->partitionLabels[0];
                defaultVals->indexToBeDeleted = 0;

                deleteNNstats();

                // Delete the vector from the front of the sliding window.
                windowValues.erase( windowValues.begin() );

                // Update the distance matrix and the NN statistics.
                updateStats();

                std::cout << "3. Point Added" << '\n';

                std::cout << "pointCounter = " << pointCounter << '\n';
                for (auto const& pair: defaultVals->avgNNDistPartitions)
                {
                    std::cout << "{" << pair.first << ": " << pair.second << "}\n";
                }
                std::cout << "==============================================================================================" << '\n';

            }
            else {   // The current vector is assigned to one of the existing partitions:
                // In this case, make sure the point to be deleted does not belong to the target partition. In particular,
                // we'll delete the oldest point from a partition label != targetPartition label.

                // Find the first occurrence of a partition label != targetPartition label.
                for (auto delIndex = 0; delIndex < defaultVals->partitionLabels.size(); delIndex++) {
                    if ( defaultVals->partitionLabels[delIndex] != defaultVals->targetPartition ) {
                        defaultVals->indexToBeDeleted = delIndex;
                        break;
                    }

                }

                defaultVals->labelToBeDeleted = defaultVals->partitionLabels[defaultVals->indexToBeDeleted];

                deleteNNstats();
                windowValues.erase( windowValues.begin() + defaultVals->indexToBeDeleted );

                // Update the distance matrix and the NN statistics.
                updateStats();

                std::cout << "4. Point Added" << '\n';

                std::cout << "pointCounter = " << pointCounter << '\n';
                for (auto const& pair: defaultVals->avgNNDistPartitions)
                {
                    std::cout << "{" << pair.first << ": " << pair.second << "}\n";
                }
                std::cout << "==============================================================================================" << '\n';
            }
        }

        return true;
    }

    return false;
}

//bool sampleStreamEvaluator(std::vector<double>& vector, std::vector<std::vector<double>>& window){
//	utils ut;
//
//	//Do some evaluation of whether the point should stay or not
//	//		For now, let's look at the deviation of connections
//
//	auto reps = ut.nearestNeighbors(vector, window);
//
//	double sum = std::accumulate(reps.begin(), reps.end(), 0.0);
//	double mean = sum / reps.size();
//
//	std::vector<double> diff(reps.size());
//	std::transform(reps.begin(), reps.end(), diff.begin(),std::bind2nd(std::minus<double>(), mean));
//	double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
//	double stdev = std::sqrt(sq_sum / reps.size());
//
//	std::sort(reps.begin(), reps.end());
//	std::vector<double> kNN;
//	int k = 20;
//
//	for(int i = 0; i < k; i++){
//		kNN.push_back(reps[i]);
//	}
//
//	double sum_NN = std::accumulate(kNN.begin(), kNN.end(), 0.0);
//	double mean_NN = sum_NN / kNN.size();
//
//	std::vector<double> diff_NN(kNN.size());
//	std::transform(kNN.begin(), kNN.end(), diff_NN.begin(),std::bind2nd(std::minus<double>(), mean_NN));
//	double sq_sum_NN = std::inner_product(diff_NN.begin(), diff_NN.end(), diff_NN.begin(), 0.0);
//	double stdev_NN = std::sqrt(sq_sum_NN / kNN.size());
//
//	if (true){//stdev_NN > 10000){
//		//std::cout << "\tAccept: (stdev > 0.5 , " << stdev << ")" << std::endl;
//		return true;
//	}
//	//std::cout << "\tReject: (stdev > 0.5 , " << stdev << ")" << std::endl;
//	return false;
//}





// runPipe -> Run the configured functions of this pipeline segment

//pipePacket slidingWindow::runPipe(pipePacket inData)
pipePacket slidingWindow::runPipe(pipePacket inData)  //basepipe.cpp
{
    utils ut;
    readInput rp;
	pPack = &inData;
	
	int cp;
	
	//Store our distance matrix
	// std::vector<std::vector<double>> distMatrix;


    // For this pipe, we construct a sub-pipeline:
    //		1. Read data vector by vector, push into slidingWindow evaluation
    //		2. IF the point is to be inserted, push into the simplexTree (LRU ordered)
    //		3. IF 100 points have passed, generate persistence intervals
    //
    // Loop this subpipeline until there's no more data

    std::vector<std::vector<double>> windowValues;

    std::vector<double> currentVector;
    if(rp.streamInit(inputFile))
    {
        // int pointCounter = 1;
        std::cout << "------------------defaukey------------------" << defaultVals->key << std::endl;

        while(rp.streamRead(currentVector))
        {
            // Initialize the sliding window. During the initialization, let's assume all points from the stream
            // belong to Partition 0.
            
             csyoutDebug("*windowValues.size*",int(windowValues.size()),true, PRCDEBUGGER);
             csyoutDebug("*pointCounter*",pointCounter,true,PRCDEBUGGER);
             csyoutDebug("*windowMaxSize",defaultVals->windowMaxSize,true, PRCDEBUGGER);
            
            if (windowValues.size() < defaultVals->windowMaxSize)
            {
                windowValues.push_back(currentVector);
                defaultVals->windowKeys.push_back(defaultVals->key);
                defaultVals->partitionLabels.push_back(defaultVals->targetPartition);
                defaultVals->key++;
                
                cp=defaultVals->key -1;
                csyoutDebug("=======csyoutDebug begin=======",cp,true, PRCDEBUGGER);
                csyoutDebug("windowKeys",defaultVals->key,true, PRCDEBUGGER);
                csyoutDebug("windowValues",windowValues[cp],PRCDEBUGGER);
                csyoutDebug("currentVector",currentVector,PRCDEBUGGER);
                csyoutDebug("partitionLabels",defaultVals->partitionLabels[cp],true,PRCDEBUGGER);
            
                
                 
                 
         
                //If we've reached window max size, generate the initial complex
                if (windowValues.size() == defaultVals->windowMaxSize)  // when read 10 data points, construct a window (672-737)
                {
                    std::cout << "Initializing complex" << std::endl;

                   inData.originalData = windowValues;
             

                    distMatrix.resize(inData.originalData.size(), std::vector<double>(inData.originalData.size(),0));
                    
                    //csyoutDebug("distMatrix",distMatrix,PRCDEBUGGER);

					//Iterate through each vector
                    for(unsigned i = 0; i < inData.originalData.size(); i++)
                    {
                        if(!inData.originalData[i].empty())
                        {
                            std::vector<double> distsFromCurrVect;

                            for(unsigned j = 0; j < inData.originalData.size(); j++) // calculate the distance between i_th point with other point (687-701)
                            {
                                if (j < i)
                                {
                                    distsFromCurrVect.push_back( distMatrix[j][i] );
                                }
                                else if (j > i)
                                {
                                    //Calculate vector distance
                                    auto dist = ut.vectors_distance(inData.originalData[i], inData.originalData[j]);
                                    distMatrix[i][j] = dist;
                                    distsFromCurrVect.push_back( dist );
                                }

                            }
                            
                            //csyoutDebug("distMatrix744 \n",distMatrix,PRCDEBUGGER);
                            
                
                            auto tempIndex = std::min_element(distsFromCurrVect.begin(), distsFromCurrVect.end()) - distsFromCurrVect.begin(); // calculate the location index of the min element in the vector
                            if (tempIndex < i)   //(704-710) find the index of the i point's nearest neighbor within the window.
                                defaultVals->nnIndices.push_back(tempIndex);   //right top triangle
                            else  // when tempIndex>=i
                                defaultVals->nnIndices.push_back(tempIndex + 1); //? why +1  left bottom triangle

                            csyoutDebug("nnIndices",defaultVals->nnIndices,PRCDEBUGGER);
                            
                            
                            csyoutDebug("i=",int(i),false, PRCDEBUGGER);
                            csyoutDebug("tempIndex=",int(tempIndex),false, PRCDEBUGGER);
                            
                          
                           
                            auto nnDistFromCurrVect = *std::min_element( distsFromCurrVect.begin(), distsFromCurrVect.end() ); // pick the nearest neighbor distance of each point within the window. (709-710)
                            defaultVals->nnDists.push_back( nnDistFromCurrVect );
                            
                            csyoutDebug("nnDists",defaultVals->nnDists,PRCDEBUGGER);
                            
                            csyoutDebug("nnIndices nnDists",defaultVals->nnIndices, defaultVals->nnDists, PRCDEBUGGER);
                        }
                    } // end of 693 for
                    
                    //distMatrix.clear(); test for simplexTree.insert() csydebug
                    inData.complex->setDistanceMatrix(&distMatrix);
                    
                    csyoutDebugComplex("inData.complex 772", inData.complex, PRCDEBUGGER);

					for(auto a : windowValues)
						inData.complex->insert();  // call insert() function from simplexTree.cpp with distMatrix(line:794 inData.complex->setDistanceMatrix(&distMatrix);)

					//csyoutDebugComplex("inData.complex 776", inData.complex, PRCDEBUGGER);

                    // Set the stream evaluator
                    inData.complex->setStreamEvaluator(&this->nnBasedEvaluator);  //this part is for setting function nnBasedEvaluator() defined in line346 which is execuated in simplexTree.cpp line 138/178
																				  // complex is in simplexBase* type; (in simplexBase.hpp) define in pipePacket.hpp
																							 
																				  //void simplexBase::setStreamEvaluator(bool (*f) (std::vector<double>&, std::vector<std::vector<double>>&))
																						 //{
																						    //streamEval = f;     
																						    //ut.writeLog(simplexType,"Changed stream evaluator");
																						    //return;
																						 //}
																								
																				 //Stream evaluator - this uses a function to determine if points should be inserted into the complex
	                                                                               //bool (*streamEval) (std::vector<double>&, std::vector<std::vector<double>>&); in simplexbase.hpp

                    std::cout << "Returning from complex initializer" << std::endl;
                    
                    //csyoutDebugComplex("inData.complex 795", inData.complex, PRCDEBUGGER);

                    // Find the average nearest neighbor distance in the existing partition (i.e. Partition 0).
                    auto avgNNDistPartition0 = std::accumulate(defaultVals->nnDists.begin(), defaultVals->nnDists.end(), 0.0) / defaultVals->nnDists.size();
                    defaultVals->avgNNDistPartitions[defaultVals->targetPartition] = avgNNDistPartition0;
                    defaultVals->numPointsPartn[defaultVals->targetPartition] = defaultVals->windowMaxSize;  // numPointsPrtn is a dictionary to store the number of points in each partition.
                    defaultVals->maxKeys[defaultVals->targetPartition] = defaultVals->key - 1;  // A dictionary to store the maxKey of each partition.  maxkey 10-1;
                    csyoutDebugComplex("inData.complex 802", inData.complex, PRCDEBUGGER);
                    
                    
                    csyoutDebug("avgNNDistPartition0 ",avgNNDistPartition0, true, PRCDEBUGGER);
					csyoutDebug("avgNNDistPartitions[defaultVals->targetPartition]", defaultVals->avgNNDistPartitions[defaultVals->targetPartition],true, PRCDEBUGGER);
					csyoutDebug("[defaultVals->targetPartition]", defaultVals->targetPartition,true,PRCDEBUGGER);
					csyoutDebug("numPointsPartn[defaultVals->targetPartition] ",defaultVals->numPointsPartn[defaultVals->targetPartition], true, PRCDEBUGGER);
					csyoutDebug("maxKeys[defaultVals->targetPartition] ", defaultVals->maxKeys[defaultVals->targetPartition], true, PRCDEBUGGER);
                    
                    
                    
                    
                    
                } //end of 672  if (windowValues.size() == defaultVals->windowMaxSize)


            } // end of  664  when if (windowValues.size() < defaultVals->windowMaxSize) is true, repeat read data until get 10 data
            
            else   // if (windowValues.size() == defaultVals->windowMaxSize)
            {
																									  //bool simplexBase::insertIterative(std::vector<double>&, std::vector<std::vector<double>>&, int&, int&){
																										//ut.writeLog(simplexType,"No insert iterative function defined");
																										//return false;
																										//}
               
              
                if(inData.complex->insertIterative(currentVector, windowValues, defaultVals->keyToBeDeleted, defaultVals->indexToBeDeleted))  // only insert one currentvector
                {
                    // inData.complex->deleteIndexRecurse( defaultVals->keyToBeDeleted );

                    // Insert the current vector, its key and partition label into the rear ends of the corresponding containers.
                    
                    csyoutDebug("*else windowValues.size*",int(windowValues.size()),true,PRCDEBUGGER);
                    csyoutDebug("*else pointCounter*",pointCounter,true,PRCDEBUGGER);
                    
                    windowValues.push_back(currentVector); 
                    defaultVals->windowKeys.push_back(defaultVals->key);
                    defaultVals->partitionLabels.push_back(defaultVals->targetPartition);
                    defaultVals->maxKeys[defaultVals->targetPartition] = defaultVals->key;  // Insert or update the maxKey of the target partition.
                    
                  
                    cp = (defaultVals->key % WINDOWSIZE);
                    
                    csyoutDebug("=======csyoutDebug else begin=======",cp,true,PRCDEBUGGER);
                    csyoutDebug("Key",defaultVals->key,true,PRCDEBUGGER);
                    
             
					//csyoutDebug("maxKeys","maxKeys",PRCDEBUGGER);
					//csyoutDebug("numPointsPartn","numPointsPartn",PRCDEBUGGER);
					//csyoutDebug("avgNNDistPartitions","avgNNDistPartitions",PRCDEBUGGER);
					
					csyoutDebug("maxKeys",defaultVals->maxKeys,PRCDEBUGGER);
					csyoutDebug("numPointsPartn",defaultVals->numPointsPartn,PRCDEBUGGER);
					csyoutDebug("avgNNDistPartitions",defaultVals->avgNNDistPartitions,PRCDEBUGGER);
				
					csyoutDebug("windowKeys",defaultVals->key,true,PRCDEBUGGER);
					csyoutDebug("windowValues",windowValues[cp],PRCDEBUGGER);
					
					csyoutDebug("windowValues 2D",windowValues,PRCDEBUGGER);
					
					csyoutDebug("currentVector",currentVector,PRCDEBUGGER);
					csyoutDebug("partitionLabels",defaultVals->partitionLabels[cp],true, PRCDEBUGGER); 
                    
                    defaultVals->key++;
                    
                    csyoutDebug("else avgNNDistPartitions[defaultVals->targetPartition]", defaultVals->avgNNDistPartitions[defaultVals->targetPartition],true, PRCDEBUGGER);
					csyoutDebug("else [defaultVals->targetPartition]", defaultVals->targetPartition,true,PRCDEBUGGER);
					csyoutDebug("else numPointsPartn[defaultVals->targetPartition] ",defaultVals->numPointsPartn[defaultVals->targetPartition], true,PRCDEBUGGER);
					csyoutDebug("else maxKeys[defaultVals->targetPartition] ", defaultVals->maxKeys[defaultVals->targetPartition],true, PRCDEBUGGER);

                }
                

                
            }  //end of 757 else 
		 

            //Check if we've gone through 100 points
            if(pointCounter % INTERVAL == 0 && pointCounter >= defaultVals->windowMaxSize) //++++++++++++++++++++++++ end of constructing a window, window size=10 ++++++++++++++++++++++++++++
            {
                // Build and trigger remaining pipeline. It should only require the computation of persistence
                // intervals from the complex being maintained.
                
                std::cout << "<<<<<<<<<<<<<<<<<<<<<<pointCounter: " << pointCounter << "\tSimplex Count: " << inData.complex->simplexCount() << "\tVertex Count: " << inData.complex->vertexCount() << std::endl;
                std::cout << "\tWindowSize: " << windowValues.size() << std::endl;
                
                inData.originalData = windowValues;
             
                
                //runSubPipeline(inData);
               std::cout <<"********************1**********pointCounter****************:"<< pointCounter-9 <<"-" << pointCounter<< "\n";
               winp++;
               
               int winpbak=winp;  
               if (winp % (WINDOWNUMS+1) == 0 ) // avoid beyond boundaries
                {			
				 winp = 1;
			
			    }

                runSubPipeline(inData,winp);  // process a window   
                
                if (winpbak % (WINDOWNUMS) == 0 ) // avoid beyond boundaries
                {
				slidingWindow::computeTotal();
				avgfirst.clear();
				avgsecond.clear();
				numPointfirst.clear();
				numPointsecond.clear();
				maxKeyfirst.clear();
				maxKeysecond.clear();
			    slidingWindow::test();
			    
			    }
               
         std::cout<<"------------------------END----------------------------------------" <<std::endl;
         std::cout<<"------------------------" <<pointCounter<<"------------------------" <<std::endl;
         std::cout<<"------------------------END----------------------------------------" <<std::endl;
        
            }
            pointCounter++;
            
            currentVector.clear();  // clear currentvector and then read next vector
            
        } //******************************** end of while(rp.streamRead(currentVector)) line 659, finish read the whole 1010 data ************
        
      
        
        //Probably want to trigger the remaining pipeline one last time... //if the whole data=1015, then this part will process the last 5 data
        if((pointCounter - 1) % INTERVAL != 0)
        {
            std::cout << "pointCounter: " << pointCounter << "\tSimplex Count: " << inData.complex->simplexCount() << "\tVertex Count: " << inData.complex->vertexCount() << std::endl;

           std::cout <<"****************2****************pointCounter********************:"<< pointCounter-9 <<"-" << pointCounter<< "\n";
            //runSubPipeline(inData);
            winp++;
            int winpbak=winp; 
             if (winp % (WINDOWNUMS+1) == 0 )   // avoid beyond boundaries
               {
				winp = 1;
			   }
			   
           runSubPipeline(inData,winp);
           
           if (winpbak % (WINDOWNUMS) == 0 )   // avoid beyond boundaries
               {
				slidingWindow::computeTotal();
				avgfirst.clear();
				avgsecond.clear();
				numPointfirst.clear();
				numPointsecond.clear();
				maxKeyfirst.clear();
				maxKeysecond.clear();
				slidingWindow::test();
			   }
           
            
        }
        
        
        
        ut.writeLog("slidingWindow", "\tSuccessfully evaluated " + std::to_string(pointCounter) + " points");
           
          writeComplexStats(inData);
     
    }
    
   

    return inData;
}

//*****************************************************************************************************************************************
//*****************************************************************************************************************************************
//*****************************************************************************************************************************************
//********************************** runPipe function finished ****************************************************************************
//*****************************************************************************************************************************************
//*****************************************************************************************************************************************
//*****************************************************************************************************************************************
//*****************************************************************************************************************************************


//void slidingWindow::writeComplexStats(pipePacket &inData)
void slidingWindow::writeComplexStats(pipePacket &inData)
{
	
    if(inData.complex->stats.size() > 30)                           //pipePacket.hpp line 17: std::string stats
    {
		
        std::ofstream file ("output/complexStats.csv");

        file << inData.complex->stats << std::endl;

        file.close();

    }
}

//void slidingWindow::runSubPipeline(pipePacket wrData)
void slidingWindow::runSubPipeline(pipePacket wrData, int winp)
{
    if(wrData.originalData.size() == 0)
        return;

    pipePacket inData = wrData;
    //outputData(inData);
   
   // outputData(inData,winp);
    
    newoutputData(inData,winp);
    
   
	std::string pipeFuncts = "rips.fast";
	// std::string pipeFuncts = "rips";
    auto lim = count(pipeFuncts.begin(), pipeFuncts.end(), '.') + 1;               // lim=1+1 when pipeFuncts = "rips.fast"
    subConfigMap["fn"] = "_" + std::to_string(repCounter);                         // slidingWindow.hpp line 16: std::map<std::string, std::string> subConfigMap;
    repCounter++;                                                                  // slidingWindow.hpp line 13: int repCounter = 0;

    //=====================generate the FastPersistence_bettis_output, comment now, will use it later=================

    //For each '.' separated pipeline function (count of '.' + 1 -> lim)
    for(unsigned i = 0; i < lim; i++)
    {
        auto curFunct = pipeFuncts.substr(0,pipeFuncts.find('.'));           // Copy to "." characters of pipeFuncts starting from position 0  = rips.
 
        pipeFuncts = pipeFuncts.substr(pipeFuncts.find('.') + 1);            // after'.' = fast

        //Build the pipe component, configure and run
        auto *cp = basePipe::newPipe(curFunct, "simplexTree");               //basePipe.cpp line 23 
                                                                             //basePipe* basePipe::newPipe(const std::string &pipeType, const std::string &complexType)

        //Check if the pipe was created and configure
        if(cp != 0 && cp->configPipe(subConfigMap))                          // slidingWindow.hpp line 50
                                                                             // bool configPipe(std::map<std::string, std::string>);
        {
            //Run the pipe function (wrapper)
            inData = cp->runPipeWrapper(inData);                             // basePipe.cpp line 51               
                                                                             // pipePacket basePipe::runPipeWrapper(pipePacket inData)
        }
        else
        {
            std::cout << "LHF : Failed to configure pipeline: " << curFunct << std::endl;
        }
    }

   
    return;
}

// configPipe -> configure the function settings of this pipeline segment
bool slidingWindow::configPipe(std::map<std::string, std::string> configMap)
{
    std::string strDebug;
    subConfigMap = configMap;

    auto pipe = configMap.find("debug");
    if(pipe != configMap.end())
    {
        debug = std::atoi(configMap["debug"].c_str());
        strDebug = configMap["debug"];
    }
    pipe = configMap.find("outputFile");
    if(pipe != configMap.end())
        outputFile = configMap["outputFile"].c_str();

    ut = utils(strDebug, outputFile);

    pipe = configMap.find("inputFile");
    if(pipe != configMap.end())
        inputFile = configMap["inputFile"].c_str();

    pipe = configMap.find("epsilon");
    if(pipe != configMap.end())
        epsilon = std::atof(configMap["epsilon"].c_str());
    else return false;

    pipe = configMap.find("dimensions");
    if(pipe != configMap.end())
    {
        dim = std::atoi(configMap["dimensions"].c_str());
    }
    else return false;

    configured = true;
    ut.writeDebug("slidingWindow","Configured with parameters { input: " + configMap["inputFile"] + ", dim: " + configMap["dimensions"] + ", eps: " + configMap["epsilon"] + ", debug: " + strDebug + ", outputFile: " + outputFile + " }");

    return true;
}

// outputData -> used for tracking each stage of the pipeline's data output without runtime
//void slidingWindow::outputData(pipePacket inData)
void slidingWindow::outputData(pipePacket inData, int winp)
{
    //std::ofstream file ("output/" + pipeType + "_" + std::to_string(repCounter) + "_output.csv");
    std::cout << "------------------------------winp--------------------------"<< winp<<"\n";
    std::ofstream file ("output/" + pipeType + "_" + std::to_string(repCounter) + "_output"+ std::to_string(winp)+".csv");

    sumpPackT[winp] = &inData;
    
    //for(auto a : inData.originalData)
    //{
        //for(auto d : a)
        //{
            //file << d << ",";
        //}
        //file << "\n";
    //}
    
       for(auto a : sumpPackT[winp]->originalData)
    {
        for(auto d : a)
        {
            file << d << ",";
        }
        file << "\n";
    }
    
    //=========================new method===============================
    memcpy(slidingWindow::sumDefaultValsT[winp], slidingWindow::defaultVals, sizeof(struct EvalParams));
    
    //=======print for defaultvals windowkeys,partitionLabels,nnindices,nnDists, avgNNDistPartitions, numPointsPartn,maxKeys, distsFromCurrVec=========
     //for(unsigned i = 0; i < defaultVals->windowMaxSize; i++)
    //{		
		//std::cout <<"******windowKeys;partitionLabels;nnIndices;nnDists;*******************"<<std::endl;	
		//std::cout << "windowKeys[" << i << "]:" << defaultVals->windowKeys[i] << std::endl;
		//std::cout << "partitionLabels["<< i <<"]:" << defaultVals->partitionLabels[i] << std::endl;
		//std::cout << "nnIndices[" << i << "]:" << defaultVals->nnIndices[i] << std::endl;
		//std::cout << "nnDists["<< i <<"]:" << defaultVals->nnDists[i] << std::endl;
			
	//}
    
    //for (auto it = defaultVals->avgNNDistPartitions.cbegin(); it!= defaultVals->avgNNDistPartitions.cend();++it)
	//{
		//std::cout <<"********************avgNNDistPartitions***************************"<<std::endl;
		//std::cout << "avgNNDistPartitions.first:" << it->first << "," << "avgNNDistPartitions.second:"<< it -> second << " " << "\n";
		
    //}
    
    
    //for (const auto& elem : defaultVals->numPointsPartn)
	//{
	  //std::cout <<"*****************************numPointsPartn**************************"<<std::endl;
	  //std::cout << "numPointsPartn.first:" << elem.first << "," << "numPointsPartn.second:"<< elem.second << " " << "\n";
	//}
	
    
    //for (unsigned j=0;j<defaultVals->distsFromCurrVec.size();j++)
	//{
		//std::cout <<"--------------------distsFromCurrVec--------------------------------"<<std::endl;
		//std::cout << "distsFromCurrVec["<< j <<"]:" << defaultVals->distsFromCurrVec[j] << std::endl;
		
     //}
    
        //for (const auto& elemk : defaultVals->maxKeys)
	//{
	  //std::cout <<"*****************************maxKeys**************************"<<std::endl;
	  //std::cout << "maxKeys.first:" << elemk.first << "," << "maxKeys.second:"<< elemk.second << " " << "\n";
	//}
    
         
    //======print for sumDefaultvals windowkeys,partitionLabels,nnindices,nnDists, avgNNDistPartitions, numPointsPartn,maxKeys, distsFromCurrVec=========
     for(unsigned i = 0; i < sumDefaultValsT[winp]->windowMaxSize; i++)
    {		
		std::cout <<"******sumDefaultVals_windowKeys;partitionLabels;nnIndices;nnDists;*******************"<<std::endl;	
		std::cout << "windowKeys[" <<"winp="<< winp <<"-i=" << i << "]:" << sumDefaultValsT[winp]->windowKeys[i] << std::endl;
		std::cout << "partitionLabels["<< i <<"]:" << sumDefaultValsT[winp]->partitionLabels[i] << std::endl;
		std::cout << "nnIndices[" << i << "]:" << sumDefaultValsT[winp]->nnIndices[i] << std::endl;
		std::cout << "nnDists["<< i <<"]:" << sumDefaultValsT[winp]->nnDists[i] << std::endl;
			
	}
    
    for (auto it = sumDefaultValsT[winp]->avgNNDistPartitions.cbegin(); it!= sumDefaultValsT[winp]->avgNNDistPartitions.cend();++it)
	{
		std::cout <<"********************sumDefaultVals_avgNNDistPartitions***************************"<<std::endl;
		std::cout << "avgNNDistPartitions.first:" << it->first << "," << "avgNNDistPartitions.second:"<< it -> second << " " << "\n";
		
    }
    
    
    for (const auto& elem : sumDefaultValsT[winp]->numPointsPartn)
	{
	  std::cout <<"*****************************sumDefaultVals_numPointsPartn**************************"<<std::endl;
	  std::cout << "numPointsPartn.first:" << elem.first << "," << "numPointsPartn.second:"<< elem.second << " " << "\n";
	}
	
	
	   for (unsigned j=0;j<sumDefaultValsT[winp]->distsFromCurrVec.size();j++)
	{
		std::cout <<"--------------------sumDefaultVals_distsFromCurrVec--------------------------------"<<std::endl;
		std::cout << "distsFromCurrVec["<< j <<"]:" << sumDefaultValsT[winp]->distsFromCurrVec[j] << std::endl;
		
     }
	
	
    
     for (auto it = sumDefaultValsT[winp]->maxKeys.cbegin(); it!= sumDefaultValsT[winp]->maxKeys.cend();++it)
	{
		std::cout <<"********************sumDefaultVals_maxkeys***************************"<<std::endl;
		std::cout << "maxKeys.first:" << it->first << "," << "maxKeys.second:"<< it -> second << " " << "\n";
		
    }
  
 
    
    file << std::endl;

    file.close();
    return;
}

//-----------------------newoutput method---------------------------
void slidingWindow::newoutputData(pipePacket inData, int wp)
{
    //std::ofstream file ("output/" + pipeType + "_" + std::to_string(repCounter) + "_output.csv");
    std::cout << "----------------------new--------winp--------------------------"<< wp<<"\n";
    std::ofstream file ("output/" + pipeType + "_" + std::to_string(repCounter) + "_newoutput"+ std::to_string(wp)+".csv");
     
    sumpPackT[wp] = &inData;
    
    sumoriginalData[wp].resize(WINDOWSIZE);
    
    for (int i=0;i< WINDOWSIZE;i++)
    {
		
		sumoriginalData[wp][i].resize(sumpPackT[wp]->originalData[0].size());
			
	}
     
    std::cout <<"sumoriginalData:" <<std::endl;
	
	 for (unsigned i = 0; i < sumoriginalData[wp].size(); i ++ ) 
	 {
        for (unsigned j = 0; j < sumoriginalData[wp][0].size(); j++ ) 
        {
            std::cout << j << ' ';
        }
        std::cout << std::endl;
    }
	std::cout << std::endl;
	std::cout << std::endl;
    
        
        file << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"<< "\n";
        int i=0;
        for(auto a : sumpPackT[wp]->originalData)
      {
		int j=0;
        for(auto d : a)
        {	
           // file << d << ",";
            sumoriginalData[wp][i][j]=d;
            file << sumoriginalData[wp][i][j] << ",";
            j++;
        }
        file << "\n";
         i++;
      }
   
       
    
    //=========================new method===============================
    memcpy(slidingWindow::sumDefaultValsT[wp], slidingWindow::defaultVals, sizeof(struct EvalParams));
    
    
    
    
         
    //======print for sumDefaultvals windowkeys,partitionLabels,nnindices,nnDists, avgNNDistPartitions, numPointsPartn,maxKeys, distsFromCurrVec=========
    
    
    
    //std::cout <<"******sumDefaultVals_windowKeys;partitionLabels;nnIndices;nnDists;*******************"<<std::endl;	
    
    if(NEWDEBUGGER==1)
  {
    
    std::cout << "windowMaxSize:" << defaultVals->windowMaxSize << std::endl;
    std::cout << "key:" << defaultVals-> key << std::endl;
    std::cout << "targetPartition: The partition membership of a new vector, if the new vector is to be added to the window." << defaultVals-> targetPartition << std::endl;
    
    std::cout << "windowKeys:" << std::endl;
    for(unsigned i = 0; i < sumDefaultValsT[wp]->windowMaxSize; i++)
    {
		std::cout <<"[" << i << "]:" << sumDefaultValsT[wp]->windowKeys[i]<< " ";
	}
    std::cout<< " " << "\n";
    std::cout << std::endl;
	
    
    std::cout << "partitionLabels:"<< std::endl;
    for(unsigned i = 0; i < sumDefaultValsT[wp]->windowMaxSize; i++)
    {
		std::cout <<"[" << i << "]:" << sumDefaultValsT[wp]->partitionLabels[i] << " ";
	}
    std::cout<< " " << "\n";
    std::cout << std::endl;
    
    std::cout << "nnIndices: A container to store the index of each point's nearest neighbor within the window."<< std::endl;
    for(unsigned i = 0; i < sumDefaultValsT[wp]->windowMaxSize; i++)
	{
		std::cout <<"[" << i << "]:" << sumDefaultValsT[wp]->nnIndices[i] << " ";
	}
    std::cout<< " " << "\n";
    std::cout << std::endl;
    
    std::cout << "nnDists: A container to store the nearest neighbor distance of each point within the window."<< std::endl;
    for(unsigned i = 0; i < sumDefaultValsT[wp]->windowMaxSize; i++)
	{
		std::cout <<"[" << i << "]:" << sumDefaultValsT[wp]->nnDists[i]<< " ";
	}
    std::cout<< " " << "\n";
    std::cout << std::endl; 
    
    
    std::cout <<"avgNNDistPartitions:"<<std::endl;
    	std::cout << "1_th:" ;
		for (auto it = sumDefaultValsT[wp]->avgNNDistPartitions.cbegin(); it!= sumDefaultValsT[wp]->avgNNDistPartitions.cend();++it)
		{
			std::cout << it->first << " ";
			avgfirst.push_back(int(it->first));
		}
		std::cout<< " " << "\n";
		
		std::cout << "2_nd:";
        for (auto it = sumDefaultValsT[wp]->avgNNDistPartitions.cbegin(); it!= sumDefaultValsT[wp]->avgNNDistPartitions.cend();++it)
		{
			std::cout << it->second << " ";	
			avgsecond.push_back(int(it->second));	
		}
        std::cout<< " " << "\n";
        std::cout << std::endl;
    
    
    std::cout <<"numPointsPartn: A dictionary to store the number of points in each partition." <<std::endl;
		std::cout << "1_th:";
		for (const auto& elem : sumDefaultValsT[wp]->numPointsPartn)
		{
			std::cout << elem.first << " ";
		}
		std::cout<< " " << "\n";
    
		std::cout << "2_nd:";
		for (const auto& elem : sumDefaultValsT[wp]->numPointsPartn)
		{
			std::cout << " " << elem.second << " ";	
		}
		std::cout<< " " << "\n";
		std::cout << std::endl;
    
       
    std::cout <<"maxKeys: A dictionary to store the maxKey of each partition."<<std::endl;
		std::cout << "1_th:" ;
		for (auto it = sumDefaultValsT[wp]->maxKeys.cbegin(); it!= sumDefaultValsT[wp]->maxKeys.cend();++it)
		{
			std::cout << it->first << " ";
		}
		std::cout<< " " << "\n";
		
		std::cout << "2_nd:";
        for (auto it = sumDefaultValsT[wp]->maxKeys.cbegin(); it!= sumDefaultValsT[wp]->maxKeys.cend();++it)
		{
			std::cout << it->second << " ";		
		}
        std::cout<< " " << "\n";
        std::cout << std::endl;
        
    
    std::cout <<"distsFromCurrVec:"<<std::endl;
		for (unsigned j=0;j<sumDefaultValsT[wp]->distsFromCurrVec.size();j++)
		{  
			std::cout <<"[" << j << "]:" << sumDefaultValsT[wp]->distsFromCurrVec[j]<< " ";
		}
		std::cout<< " " << "\n";
		std::cout << std::endl;
     
	  std::cout << "keyToBeDeleted:" << defaultVals->keyToBeDeleted << std::endl;
	  std::cout << "labelToBeDeleted:" << defaultVals->labelToBeDeleted << std::endl;
	  std::cout << "indexToBeDeleted:" << defaultVals->indexToBeDeleted << std::endl;
	  std::cout << "nnDistToBeDeleted:" << defaultVals->nnDistToBeDeleted << std::endl;
	
	
	//==============print for ppack===========
	
	
	std::cout << "originalLabels: " ;
	std::cout << "size =" << sumpPackT[wp]->originalLabels.size()<< std::endl;    // why is 0?

	
	for (unsigned i=0; i=sumpPackT[wp]->originalLabels.size();i++)
	{
		std::cout <<"[" << i << "]:" << sumpPackT[wp]->originalLabels[i] << " ";
	}
    std::cout<< " " << "\n";  
    std::cout << std::endl;
	
	std::cout << "fulldata:" ; 
	std::cout << "size =" << sumpPackT[wp]->fullData.size()<< std::endl; 
	 for (unsigned i = 0; i < sumpPackT[wp]->fullData.size(); i ++ ) 
	 {
        for (unsigned j = 0; j < sumpPackT[wp]->fullData[0].size(); j++ ) 
        {
            std::cout << sumpPackT[wp]->fullData[i][j] << ' ';
        }
        std::cout << std::endl;
    }
	std::cout << std::endl;
	std::cout << std::endl;
	
	
	std::cout <<"originalData:" <<std::endl;
	
	 for (unsigned i = 0; i < sumpPackT[wp]->originalData.size(); i ++ ) 
	 {
        for (unsigned j = 0; j < sumpPackT[wp]->originalData[0].size(); j++ ) 
        {
            std::cout << sumpPackT[wp]->originalData[i][j] << ' ';
        }
        std::cout << std::endl;
    }
	std::cout << std::endl;
	std::cout << std::endl;
	
	
	
	std::cout <<"distMatrix:" <<std::endl;
	std::cout << "size =" << sumpPackT[wp]->distMatrix.size()<< std::endl; 
	 for (unsigned i = 0; i < sumpPackT[wp]->distMatrix.size(); i ++ ) 
	 {
        for (unsigned j = 0; j < sumpPackT[wp]->distMatrix[0].size(); j++ ) 
        {
            std::cout << sumpPackT[wp]->distMatrix[i][j] << ' ';
        }
        std::cout << std::endl;
    }
	std::cout << std::endl;
	std::cout << std::endl;
	
	
	std::cout <<"boundaries:" <<std::endl;
	std::cout << "size =" << sumpPackT[wp]->boundaries.size()<< std::endl; 
	 for (unsigned i = 0; i < sumpPackT[wp]->boundaries.size(); i ++ ) 
	 {
        for (auto d: sumpPackT[wp]->boundaries[i]) 
        {
            std::cout << d << ' ';
        }
        std::cout << std::endl;
    }
	std::cout << std::endl;
	std::cout << std::endl;
	
	
	std::cout << "bettiOutput:" << sumpPackT[wp]->bettiOutput << std::endl;
	std::cout << "stats:" << sumpPackT[wp]->stats << std::endl;
	 
	
	
	std::cout <<"weights:"<<std::endl;
	
	for (auto d: sumpPackT[wp]->weights) 
        {
            std::cout << d << ' ';
        }
        std::cout << std::endl;
    
	
	
    std::cout << "complex nodeCount: "  << sumpPackT[wp]->complex->nodeCount <<" Total number of nodes stored" << std:: endl;
	std::cout << "indexCounter: "  << sumpPackT[wp]->complex->indexCounter <<" Current insertion index" << std:: endl;
	std::cout << "root: "  << sumpPackT[wp]->complex->root << " Root of the simplexNode tree (if applicable)"<< std:: endl;
	std::cout << "head: "  << sumpPackT[wp]->complex->head <<  " Root of the simplexNode tree (if applicable)" << std:: endl;
	std::cout << "maxEpsilon: "  << sumpPackT[wp]->complex->maxEpsilon <<" Maximum epsilon, loaded from configuration"<<std:: endl;
	std::cout << "maxDimension: "  << sumpPackT[wp]->complex->maxDimension << " Maximum dimension, loaded from configuration" <<std:: endl;
  }


	
	//=========================================when newdebugger=0===========================================
	
	//std::cout <<"******sumDefaultVals_windowKeys;partitionLabels;nnIndices;nnDists;*******************"<<std::endl;	
    
    if(NEWDEBUGGER==0)
  {
    
    file << "windowMaxSize:" << defaultVals->windowMaxSize << "\n";
    file<< "key:" << defaultVals-> key << "\n";
    file << "targetPartition: " << defaultVals-> targetPartition <<"(The partition membership of a new vector, if the new vector is to be added to the window.)"<< "\n";    
      
    file << "windowKeys:" << "\n";
    for(unsigned i = 0; i < sumDefaultValsT[wp]->windowMaxSize; i++)
    {
		file <<"[" << i << "]:" << sumDefaultValsT[wp]->windowKeys[i]<< " ";
	}
    file<< " " << "\n";
    file << std::endl;
	
    
    file << "partitionLabels:"<< "\n";
    for(unsigned i = 0; i < sumDefaultValsT[wp]->windowMaxSize; i++)
    {
		file <<"[" << i << "]:" << sumDefaultValsT[wp]->partitionLabels[i] << " ";
	}
    file<< " " << "\n";
    file << "\n";
    
    file << "nnIndices: (A container to store the index of each point's nearest neighbor within the window.)"<< "\n";
    for(unsigned i = 0; i < sumDefaultValsT[wp]->windowMaxSize; i++)
	{
		file <<"[" << i << "]:" << sumDefaultValsT[wp]->nnIndices[i] << " ";
	}
    file<< " " << "\n";
    file << "\n";
    
    file << "nnDists: (A container to store the nearest neighbor distance of each point within the window.)"<< "\n";
    for(unsigned i = 0; i < sumDefaultValsT[wp]->windowMaxSize; i++)
	{
		file <<"[" << i << "]:" << sumDefaultValsT[wp]->nnDists[i]<< " ";
	}
    file<< " " << "\n";
    file << "\n"; 
    
    
    file <<"avgNNDistPartitions:(the average nearest neighbor distance in the existing partition)"<<"\n";
    	file << "1_th:" ;
		for (auto it = sumDefaultValsT[wp]->avgNNDistPartitions.cbegin(); it!= sumDefaultValsT[wp]->avgNNDistPartitions.cend();++it)
		{
			file << it->first << " ";
			avgfirst.push_back(int(it->first));
		}
		file<< " " << "\n";
		
		file << "2_nd:";
        for (auto it = sumDefaultValsT[wp]->avgNNDistPartitions.cbegin(); it!= sumDefaultValsT[wp]->avgNNDistPartitions.cend();++it)
		{
			file << it->second << " ";	
			avgsecond.push_back(double(it->second));	
		}
        file<< " " << "\n";
        file << "\n";
    

    
    file <<"numPointsPartn: (A dictionary to store the number of points in each partition.)" <<"\n";
		file << "1_th:";
		for (const auto& elem : sumDefaultValsT[wp]->numPointsPartn)
		{
			file << elem.first << " ";
			numPointfirst.push_back(int(elem.first));
		}
		file<< " " << "\n";
    
		file << "2_nd:";
		for (const auto& elem : sumDefaultValsT[wp]->numPointsPartn)
		{
			file << " " << elem.second << " ";	
			numPointsecond.push_back(int(elem.second));
		}
		file<< " " << "\n";
		file << "\n";
    
       
    file <<"maxKeys: (A dictionary to store the maxKey of each partition.)"<<"\n";
		file << "1_th:" ;
		for (auto it = sumDefaultValsT[wp]->maxKeys.cbegin(); it!= sumDefaultValsT[wp]->maxKeys.cend();++it)
		{
			file << it->first << " ";
			maxKeyfirst.push_back(it->first);
		}
		file<< " " << "\n";
		
		file << "2_nd:";
        for (auto it = sumDefaultValsT[wp]->maxKeys.cbegin(); it!= sumDefaultValsT[wp]->maxKeys.cend();++it)
		{
			file << it->second << " ";		
			maxKeysecond.push_back(it->second);
		}
        file<< " " << "\n";
        file << "\n";
        
    
    file <<"distsFromCurrVec: (Compute the distances from the current vector to the existing ones in the window)"<<"\n";
		for (unsigned j=0;j<sumDefaultValsT[wp]->distsFromCurrVec.size();j++)
		{  
			file <<"[" << j << "]:" << sumDefaultValsT[wp]->distsFromCurrVec[j]<< " ";
		}
		file<< " " << "\n";
		file << "\n";
     
	  file << "keyToBeDeleted:" << defaultVals->keyToBeDeleted << "\n";
	  file << "labelToBeDeleted:" << defaultVals->labelToBeDeleted << "\n";
	  file << "indexToBeDeleted:" << defaultVals->indexToBeDeleted << "\n";
	  file << "nnDistToBeDeleted:" << defaultVals->nnDistToBeDeleted << "\n";
	
	
	//==============print for ppack===========
	
	
	file << "originalLabels: " ;
	file << "size =" << sumpPackT[wp]->originalLabels.size()<< "\n";    // why is 0?

	
	for (unsigned i=0; i=sumpPackT[wp]->originalLabels.size();i++)
	{
		file <<"[" << i << "]:" << sumpPackT[wp]->originalLabels[i] << " ";
	}
    file<< " " << "\n";  
    file << "\n";
	
	file << "fulldata:" ; 
	file << "size =" << sumpPackT[wp]->fullData.size()<< "\n"; 
	 for (unsigned i = 0; i < sumpPackT[wp]->fullData.size(); i ++ ) 
	 {
        for (unsigned j = 0; j < sumpPackT[wp]->fullData[0].size(); j++ ) 
        {
            file << sumpPackT[wp]->fullData[i][j] << ' ';
        }
        file << "\n";
    }
	file << "\n";
	file << "\n";
	
	
	file <<"originalData:" <<"\n";
	
	 for (unsigned i = 0; i < sumpPackT[wp]->originalData.size(); i ++ ) 
	 {
        for (unsigned j = 0; j < sumpPackT[wp]->originalData[0].size(); j++ ) 
        {
            file << sumpPackT[wp]->originalData[i][j] << ' ';
        }
        file << "\n";
    }
	file << "\n";
	file << "\n";
	
	
	
	file <<"distMatrix:" <<"\n";
	file << "size =" << sumpPackT[wp]->distMatrix.size()<< "\n"; 
	 for (unsigned i = 0; i < sumpPackT[wp]->distMatrix.size(); i ++ ) 
	 {
        for (unsigned j = 0; j < sumpPackT[wp]->distMatrix[0].size(); j++ ) 
        {
            file << sumpPackT[wp]->distMatrix[i][j] << ' ';
        }
        file << "\n";
    }
	file << "\n";
	file << "\n";
	
	
	file <<"boundaries:" <<"\n";
	file << "size =" << sumpPackT[wp]->boundaries.size()<< "\n"; 
	 for (unsigned i = 0; i < sumpPackT[wp]->boundaries.size(); i ++ ) 
	 {
        for (auto d: sumpPackT[wp]->boundaries[i]) 
        {
            file << d << ' ';
        }
        file << "\n";
    }
	file << "\n";
	file << "\n";
	
	
	file << "bettiOutput:" << sumpPackT[wp]->bettiOutput << "\n";
	file << "stats:" << sumpPackT[wp]->stats << "\n";
	 
	
	
	file <<"weights:"<<"\n";
	
	for (auto d: sumpPackT[wp]->weights) 
        {
            file << d << ' ';
        }
        file << "\n";
    
	
	
    file << "complex nodeCount: "  << sumpPackT[wp]->complex->nodeCount <<" (Total number of nodes stored)" << "\n";
	file << "indexCounter: "  << sumpPackT[wp]->complex->indexCounter <<" (Current insertion index)" << "\n";
	file << "root: "  << sumpPackT[wp]->complex->root << " (Root of the simplexNode tree (if applicable))"<< "\n";
	file << "head: "  << sumpPackT[wp]->complex->head <<  " (Root of the simplexNode tree (if applicable))" << "\n";
	file << "maxEpsilon: "  << sumpPackT[wp]->complex->maxEpsilon <<" (Maximum epsilon, loaded from configuration)"<< "\n";
	file << "maxDimension: "  << sumpPackT[wp]->complex->maxDimension << " (Maximum dimension, loaded from configuration)" << "\n";
}	

   
 	
	
    file << "\n";

    file.close();
    return;
}

void slidingWindow::test()
{
	std::vector<double> currentVector;
	int winp1=1;
	int winp2=2;
	int winp3=3;
	int winp4=4;
	
	//std::vector<double> currentVector=sumpPackT[winp1]->originalData[1];
	//memcpy(&currentVector, &sumpPackT[winp1]->originalData[1],sumpPackT[winp1]->originalData[0].size());
	
	//std::vector<double> distsFromCurrVec = ut.nearestNeighbors(sumoriginalData[winp1][1], sumoriginalData[winp2]);
	std::vector<double> distsFromCurrVec1 = ut.nearestNeighbors(sumoriginalData[winp1][1], sumoriginalData[winp1]);
	
	
	//std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~distsFromCurrVec:"<< "\n";
    //for(unsigned i = 0; i < distsFromCurrVec.size(); i++)
    //{
		//std::cout <<"[" << i << "]:" << distsFromCurrVec[i] << " ";
	//}
    //std::cout << " " << "\n";
    //std::cout << "\n";
    
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~1~~~~~~~~~~~~~~~~~~~~~~~~~~~~~distsFromCurrVec1:"<< "\n";
    for(unsigned i = 0; i < distsFromCurrVec1.size(); i++)
    {
		std::cout <<"[" << i << "]:" << distsFromCurrVec1[i] << " ";
	}
    std::cout << " " << "\n";
    std::cout << "\n";
    
    
    
    
    
    //===========================================test computing==========================  
    int k1=0;
    
    totalDefaultVals->key = WINDOWNUMS*WINDOWSIZE;
	totalDefaultVals->targetPartition = 0;
    
	std::cout<< "Sosososososososk1:" << winp4 << "windowMaxSize:" << sumDefaultValsT[1]->windowMaxSize<<" ";
	std::cout << "sumDefaultValsT[" << winp4 << "]" << "windowKeys[" << "0" <<"]" << sumDefaultValsT[ winp4]->windowKeys[0] << " ";
	std::cout <<"sumDefaultValsT["<<  winp4 <<"]->partitionLabels["<< 0 <<"]" << sumDefaultValsT[ winp4]->partitionLabels[0] <<std::endl;
	
	std::cout << "sumDefaultValsT[" << winp4 << "]" << "windowKeys[" <<0 <<"]" << sumDefaultValsT[winp4]->windowKeys[0] << " ";
	
	std::cout << "totalDefaultVals->windowKeys: " << " ";
	
	for (int winptr=1;winptr<=WINDOWNUMS;winptr++)           //sumDefault start from 1 to 4
	{
		for(int j = 0; j < sumDefaultValsT[winptr]->windowMaxSize; j++)
		{
			
			totalDefaultVals->windowKeys.push_back(sumDefaultValsT[winptr]->windowKeys[j]);
			//std::cout << "sumDefaultValsT[" << winptr << "]" << "windowKeys[" <<j <<"]" << sumDefaultValsT[winptr]->windowKeys[j] << " ";
			std::cout << "[" << k1 <<"]" << totalDefaultVals->windowKeys[k1] << " ";
			k1++;
		
		}
		
	
	}
	
	
	//=======================test for setXOR function in utils, remember to use "utils.functionname" to call. ut is instance of utils
    //std::set<unsigned> setA={2,4,5,6,7,8,9,10};
	//std::set<unsigned> setB={5,6,1,12,24,9};
	
	//std::vector<unsigned> v1{2,4,5,6,7,8,9,10};
	//std::vector<unsigned> v2{5,6,1,12,24,9};
	
	//auto a= ut.setXOR(setA, setB);
	 //auto b=ut.setIntersect(setA, setB, true);	
	//std::cout <<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<"\n";
	//for(auto z : c){   //print set unsigned a
			//std::cout << z << ",";
	//}
	//std::cout << "\n";
	//return;
	
	
	 
	

	
	
	
	//std::cout <<  "test xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxdistsFromCurrVec="  << distsFromCurrVec <<std::endl;
	
	//std::cout <<  "test xxxxxxxxxxxxxxxxxxxxxxxxxsumpPackT[winp1]->originalData[1]= " << std::endl;
		
	
}


void slidingWindow::csyoutDebug(std::string s, double d, bool b, bool deg)
{
	if (deg==false) return;
	if(b==true)
	{
	 std::cout<< s << "=" << d <<"\n";	
	}
	else
	{
	std::cout<< s << "=" << d <<" ";	
	}
	
}


void slidingWindow::csyoutDebug(std::string s, int d, bool b , bool deg)
{
	if (deg==false) return;
	if(b==true)
	{
		
	 std::cout<< s << "=" << d <<"\n";	
	}
	else
	{
	std::cout<< s << "=" << d <<" ";	
	}
	
}

void slidingWindow::csyoutDebug(std::string s, std::vector<double> v, bool deg)
{
	if (deg==false) return;
	std::cout<< s << ":" ;

	for (unsigned i=0;i< v.size();i++)
	{	
	  std::cout << v[i] <<" ";	
	}
	std::cout<<" " <<std::endl;

}


void slidingWindow::csyoutDebug(std::string s, std::vector<int> v, bool deg)
{
	if (deg==false) return;
	std::cout<< s << ":" ;

	for (unsigned i=0;i< v.size();i++)
	{	
	  std::cout << v[i] <<" ";	
	}
	std::cout<<" " <<std::endl;

}



void slidingWindow::csyoutDebug(std::string s, std::vector<int> v, std::vector<double> v2, bool deg)
{
	if (deg==false) return;
	std::cout<< s << ":" ;

	for (unsigned i=0;i< v.size();i++)
	{	
	  std::cout << "[ "<< i<< ", " << v[i] << " ]=" <<v2[i] <<" ";	
	}
	std::cout<<" " <<std::endl;

}







//void slidingWindow::csyoutDebug(std::string s, std::string m, bool b)
//{
	//std::cout<< s << ":" ;
	
	
	////if(compare(m,"avgNNDistPartitions")==0)
	//if((m.compare("avgNNDistPartitions"))==0)
	
	//{

		   //for (auto it = defaultVals->avgNNDistPartitions.cbegin(); it!= defaultVals->avgNNDistPartitions.cend();++it)
		//{
			//std::cout << "avgNNDistPartitions.first:" << it->first << "," << "avgNNDistPartitions.second:"<< it -> second << " " << "\n";
			
		//}
    //}
    
    ////if(compare(m,"numPointsPartn")==0)
    //if((m.compare("numPointsPartn"))==0)
    //{
		//for (const auto& elem : defaultVals->numPointsPartn)
		//{
		  //std::cout << "numPointsPartn.first:" << elem.first << "," << "numPointsPartn.second:"<< elem.second << " " << "\n";
		//}
    //}
	
	
	////if(compare(m,"maxKeys")==0)
	//if((m.compare("maxKeys"))==0)
   //{ 
		 //for (auto it = defaultVals->maxKeys.cbegin(); it!= defaultVals->maxKeys.cend();++it)
		//{
			//std::cout << "maxKeys.first:" << it->first << "," << "maxKeys.second:"<< it -> second << " " << "\n";
			
		//}
	
    //}
	
	
//}


void slidingWindow::csyoutDebug(std::string s, std::unordered_map<int, int> m, bool deg)
{
	if (deg==false) return;
	std::cout<< s << ":" ;
	
			   for (auto it = m.cbegin(); it!= m.cend();++it)
			{
				std::cout << "first:" << it->first << "," << "second:"<< it -> second << " " ;
				
			}
   		std::cout<< " "<< std::endl;
}

void slidingWindow::csyoutDebug(std::string s, std::unordered_map<int, double> m, bool deg)
{
	if (deg==false) return;
	std::cout<< s << ":" ;
	
			   //for (auto it = m.cbegin(); it!= m.cend();++it)
			//{
				//std::cout << "first:" << it->first << "," << "second:"<< it -> second << " " << "\n";
				
			//}
   
          for (const auto& elem : m)
		{
		  std::cout << "first:" << elem.first << "," << "second:"<< elem.second << " " ;
		}
		std::cout<< " "<< std::endl;
   
}


void slidingWindow::csyoutDebug(std::string s,  std::vector<std::vector<double>> v, bool deg)
{
	if (deg==false) return;
	std::cout<< s << ":" ;
	
	 for (unsigned i = 0; i < v.size(); i ++ ) 
	 {
        for (unsigned j = 0; j < v[0].size(); j++ ) 
        {
            std::cout << v[i][j] << ' ';
        }
        std::cout << std::endl;
    }
	std::cout << std::endl;
	std::cout << std::endl;	
	
}


void slidingWindow::csyoutDebugComplex(std::string s,  simplexBase* newcomplex, bool deg)
{
     if (deg==false) return;
	
	std::cout<< s << ":" ;
	std::cout<< "nodeCount:" << newcomplex->nodeCount << "\t";
	std::cout <<"simplexOffset:"<<newcomplex->simplexOffset << "\t";
	std::cout <<"indexCounter:" <<newcomplex->indexCounter << "\t";
	std::cout<< "maxEpsilon:" <<newcomplex->maxEpsilon << "\t";
	std::cout << "maxDimension:" <<newcomplex->maxDimension << "\n";
	
	
	std::cout<< "distMatrix: \n" ;
	
	//while (!newcomplex->distMatrix->empty())
	//{
		 //for (std::vector<unsigned>::iterator p = newcomplex->distMatrix->front().begin();p!= newcomplex->distMatrix->front().end();p++ )
			//{
			 //std::cout <<*p <<" ";
			//}
		 //std::cout << std::endl;
		 //newcomplex->distMatrix->erase(newcomplex->distMatrix().begin());
 	 
	 //}
			

	 for (double i = 0; i < newcomplex->distMatrix -> size(); i ++) 
	 {
        for (double j = 0; j < newcomplex->distMatrix->operator[](i).size(); j++) 
        {
            //std::cout << *(newcomplex->distMatrix)[i][j] << ' ';
            std::cout << newcomplex->distMatrix->operator []( i )[j] << ' ';
        }
        std::cout << std::endl;
    }
	std::cout << std::endl;
	std::cout << std::endl;	
	
	
	std::cout << "runningVectorIndices;";
   	for (unsigned i=0;i< newcomplex->runningVectorIndices.size();i++)
	{	
	  std::cout << "[ "<< i<< "] " << newcomplex->runningVectorIndices[i] <<" ";	
	}
	std::cout<<" " <<std::endl;
	
	std::cout <<"runningVectorCount:" << newcomplex->runningVectorCount;
	std::cout<<" " <<std::endl;
	
}

void slidingWindow::computeTotal()
{
	
	//1. compute the total default 
	int k1=0;
	int k2=0;
	int k3=0;
	int k4=0;
	int k5=0;
	
    totalDefaultVals->key = WINDOWNUMS*WINDOWSIZE;
	totalDefaultVals->targetPartition = 0;
    
	std::cout << "totalDefaultVals->windowKeys: " << "\n ";
	for (int winptr=1;winptr<=WINDOWNUMS;winptr++)           //sumDefault start from 1 to 4
	{
	    
	    totalDefaultVals->windowMaxSize = totalDefaultVals->windowMaxSize + sumDefaultValsT[winptr]->windowMaxSize;
		for(int j = 0; j < sumDefaultValsT[winptr]->windowMaxSize; j++)
		{
			
			totalDefaultVals->windowKeys.push_back(sumDefaultValsT[winptr]->windowKeys[j]);
			std::cout << "[" << k1 <<"] " << totalDefaultVals->windowKeys[k1] << " ";
			k1++;
		}
		
	}
	std::cout << "\n" ;
	
	std::cout << "totalDefaultVals->partitionLabels: " << "\n";
	for (int winptr=1;winptr<=WINDOWNUMS;winptr++)   
	 {
		 
	    for(int j = 0; j < sumDefaultValsT[winptr]->windowMaxSize; j++)
		{
		  totalDefaultVals->partitionLabels.push_back(sumDefaultValsT[winptr]->partitionLabels[j]);
		  std::cout << "[" << k2 << "] " <<  totalDefaultVals->partitionLabels[k2] << " ";
		  k2++;	  
		}
     }
	std::cout << "\n" ;
	
	std::cout << "totalDefaultVals->nnIndices: " << "\n";
	for (int winptr=1;winptr<=WINDOWNUMS;winptr++)  
	{
		for(int j = 0; j < sumDefaultValsT[winptr]->windowMaxSize; j++)
		{
			totalDefaultVals->nnIndices.push_back(sumDefaultValsT[winptr]->nnIndices[j]);
			std::cout << "[" << k3 << "] " <<  totalDefaultVals->nnIndices[k3] << " ";
			k3++;
		}
	}
	std::cout << "\n" ;
	
	std::cout << "totalDefaultVals->nnDists: " << "\n";
	for (int winptr=1;winptr<=WINDOWNUMS;winptr++)  
	{
		for(int j = 0; j < sumDefaultValsT[winptr]->windowMaxSize; j++)
		{
			totalDefaultVals->nnDists.push_back(sumDefaultValsT[winptr]->nnDists[j]);
			std::cout << "[" << k4 << "] " <<  totalDefaultVals->nnIndices[k4] << " ";
			k4++;
		}
  
	}
	std::cout << "\n" ;
	
	std::cout << "totalDefaultVals->distsFromCurrVec" << "\n";
	for (int winptr=1;winptr<=WINDOWNUMS;winptr++)
	{
	    for (unsigned j=0;j<sumDefaultValsT[winptr]->distsFromCurrVec.size();j++)
		{  
			  totalDefaultVals->distsFromCurrVec.push_back(sumDefaultValsT[winptr]->distsFromCurrVec[j]);
			  std::cout <<"[" << k5 << "] " << totalDefaultVals->distsFromCurrVec[k5] << " ";
			  k5++;
		}	
	}
	std::cout << "\n" ;
	
	
	
	
	std::cout <<"????????????????????????????????"<<" \n";
	
	std::cout << "avgsize:" << avgfirst.size() << ","<< avgsecond.size()<< "\n"; 
    for (int i=0;i<avgfirst.size();i++) 
	{
		
	std::cout << "avg["<< i<< "]:"<< avgfirst[i] << "," <<  avgsecond[i]<<"\n";	
	totalDefaultVals->avgNNDistPartitions.insert({avgfirst[i],avgsecond[i]});
	//totalDefaultVals->avgNNDistPartitions.insert(std::make_pair<int,double> (int(avgfirst[i]),double(avgsecond[i])));
    }
        std::cout <<"totalDefaultVals->avgNNDistPartitions:"<<std::endl;
    	std::cout << "1_th:" ;
		for (auto it = totalDefaultVals->avgNNDistPartitions.cbegin(); it!= totalDefaultVals->avgNNDistPartitions.cend();++it)
		{
			std::cout << it->first << " ";
		}
		std::cout<< " " << "\n";
		
		std::cout << "2_nd:";
        for (auto it = totalDefaultVals->avgNNDistPartitions.cbegin(); it!= totalDefaultVals->avgNNDistPartitions.cend();++it)
		{
			std::cout << it->second << " ";		
		}
        std::cout<< " " << "\n";
        std::cout << std::endl;
        
        totalDefaultVals->avgNNDistPartitions.clear();                  
        
        
        
    std::cout << "numPointsize:" << numPointfirst.size() << "," << numPointsecond.size()<< "\n";   
    for (int i=0;i<numPointfirst.size();i++) 
	{
		
	std::cout << "numPoint["<< i<< "]:"<< numPointfirst[i] << "," << numPointsecond[i]<<"\n";	
	totalDefaultVals->numPointsPartn.insert({numPointfirst[i],numPointsecond[i]});
    }  
     std::cout <<"totalDefaultVals->numPointsPartn:"<<std::endl;
    	std::cout << "1_th:" ;
		for (auto it = totalDefaultVals->numPointsPartn.cbegin(); it!= totalDefaultVals->numPointsPartn.cend();++it)
		{
			std::cout << it->first << " ";
		}
		std::cout<< " " << "\n";
		
		std::cout << "2_nd:";
        for (auto it = totalDefaultVals->numPointsPartn.cbegin(); it!= totalDefaultVals->numPointsPartn.cend();++it)
		{
			std::cout << it->second << " ";		
		}
        std::cout<< " " << "\n";
        std::cout << std::endl;
        
        totalDefaultVals->numPointsPartn.clear();
    
    
    
    for (int i=0;i<maxKeyfirst.size();i++) 
	{
		
	std::cout << "maxKeys["<< i<< "]:"<< maxKeyfirst[i] << "," << maxKeysecond[i]<<"\n";	
	totalDefaultVals->maxKeys.insert({maxKeyfirst[i],maxKeysecond[i]});
    }  
     std::cout <<"totalDefaultVals->maxKeys:"<<std::endl;
    	std::cout << "1_th:" ;
		for (auto it = totalDefaultVals->maxKeys.cbegin(); it!= totalDefaultVals->maxKeys.cend();++it)
		{
			std::cout << it->first << " ";
		}
		std::cout<< " " << "\n";
		
		std::cout << "2_nd:";
        for (auto it = totalDefaultVals->maxKeys.cbegin(); it!= totalDefaultVals->maxKeys.cend();++it)
		{
			std::cout << it->second << " ";		
		}
        std::cout<< " " << "\n";
        std::cout << std::endl;
        
        totalDefaultVals->maxKeys.clear();
        
  //2. compute sumorignialdata
  
  std::cout <<"originalData:" <<std::endl;
  std::cout<<"sumoriginalData[1].size: " << sumoriginalData[1].size() << " ";
  std::cout<<"sumoriginalData[1][0].size: " << sumoriginalData[1][0].size() << " " <<std::endl;
  
	
	 for (unsigned wpp = 1; wpp <= WINDOWNUMS; wpp++ ) 
	 {
        for (unsigned i = 0; i < sumoriginalData[wpp].size(); i++ )   // windowsize=10
        {
            
          for (unsigned j = 0; j < sumoriginalData[wpp][0].size(); j++ )   
           {  
            std::cout << sumoriginalData[wpp][i][j] << ",";
		   }
		  std::cout << std::endl;
        }
        std::cout << std::endl;
    }
	std::cout << std::endl;     
        
 
        
    
	////totalDefaultVals->keyToBeDeleted= defaultVals->keyToBeDeleted;
	////totalDefaultVals->labelToBeDeleted= defaultVals->labelToBeDeleted;
	////totalDefaultVals->indexToBeDeleted= defaultVals->indexToBeDeleted;
	////totalDefaultVals->nnDistToBeDeleted=defaultVals->nnDistToBeDeleted;
		
}



