/*
 * streamVR hpp + cpp extend the basePipe class for
 *
 */

#include <string>
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


int pointCounter = 1;
std::vector<std::vector<double>> distMatrix;

// basePipe constructor
template<typename nodeType>
slidingWindow<nodeType>::slidingWindow()
{
    windowMaxSize = 100;
    key = 0;
    targetPartition = 0;
    pipeType = "SlidingWindow";
    return;
}

template<typename nodeType>
void slidingWindow<nodeType>::deleteNNstats()
{
    keyToBeDeleted = windowKeys[indexToBeDeleted];

    windowKeys.erase( windowKeys.begin() + indexToBeDeleted );
    partitionLabels.erase( partitionLabels.begin() + indexToBeDeleted );

    nnIndices.erase( nnIndices.begin() + indexToBeDeleted );

    nnDistToBeDeleted = nnDists[indexToBeDeleted];
    nnDists.erase( nnDists.begin() + indexToBeDeleted );

    distsFromCurrVec.erase( distsFromCurrVec.begin() + indexToBeDeleted );

    numPointsPartn[labelToBeDeleted] = numPointsPartn[labelToBeDeleted] - 1;

    // If there are no more points left in the partition from which the deletion took place:
    if (numPointsPartn[labelToBeDeleted] == 0)
    {
        avgNNDistPartitions.erase(labelToBeDeleted);
        numPointsPartn.erase(labelToBeDeleted);
        maxKeys.erase(labelToBeDeleted);
    }

    return;

}

template<typename nodeType>
void slidingWindow<nodeType>::updateStats()
{
    // Delete the corresponding row from the distance matrix.
    // pPack->complex->distMatrix->erase(pPack->complex->distMatrix->begin() + indexToBeDeleted);
    distMatrix.erase(distMatrix.begin() + indexToBeDeleted);


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
        // pPack->complex->distMatrix->at(i).erase( pPack->complex->distMatrix->at(i).begin() + indexToBeDeleted );
        distMatrix[i].erase( distMatrix[i].begin() + indexToBeDeleted );

        // Add the new distance value to the end of each row.
        // pPack->complex->distMatrix->at(i).push_back( distsFromCurrVec[i] );
        distMatrix[i].push_back( distsFromCurrVec[i] );

        // Update NN statistics for only those partitions from which the point was deleted or to which the new point is to be added.
        // Case 1: The i-th point belongs to the partition the last point was deleted from, but not to the partition the new point is
        // to be added to. And, the point that was deleted was the nearest neighbor of the i-th point.
        if ( partitionLabels[i] == labelToBeDeleted
                &&  partitionLabels[i] != targetPartition
                && nnIndices[i] == indexToBeDeleted )
        {

            // If only one point remains in the partition from which the last point was deleted:
            if ( numPointsPartn[labelToBeDeleted] == 1 )
            {
                avgNNDistPartitions[labelToBeDeleted] = -1;
                nnIndices[i] = -1;
                nnDists[i] = -1;

            }
            else
            {
                sumDelOldNNdists = sumDelOldNNdists + nnDists[i];

                std::vector<double> memberDistsFromVect;
                std::vector<int> memberIndices;

                for(unsigned int j = 0; j < distMatrix.size(); j++)
                {
                    if ( partitionLabels[j] == labelToBeDeleted )
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
                nnIndices[i] = memberIndices[tempIdx];

                auto newNNdistFromVect = *std::min_element( memberDistsFromVect.begin(), memberDistsFromVect.end() );
                nnDists[i] = newNNdistFromVect;
                sumDelNewdists = sumDelNewdists + newNNdistFromVect;
            }

        }

        // Case 2: The i-th point belongs to the partition the new point is to be added to, but not to the partition the last point
        // was deleted from. In this case, is it possible that the point that was deleted was the nearest neighbor of the i-th point?
        else if ( partitionLabels[i] != labelToBeDeleted && partitionLabels[i] == targetPartition )
        {

            // If, previously, there was only one point in the target partition:
            if ( numPointsPartn[targetPartition] == 1 )
            {

                // Update the NN statistics for the "new" two-point partition.
                nnIndices[i] = windowMaxSize - 1;
                nnIndices.push_back(i);

                nnDists[i] = distsFromCurrVec[i];
                nnDists.push_back(nnDists[i]);

                avgNNDistPartitions[targetPartition] = nnDists[i];
                numPointsPartn[targetPartition] = 2;
                twoPointPartition = 1;
            }

            else
            {

                distsFromCurrVecTP.push_back( distsFromCurrVec[i] );
                tpIndices.push_back(i);

                if ( distsFromCurrVec[i] < nnDists[i] )
                {
                    sumAddOldNNdists = sumAddOldNNdists + nnDists[i];

                    nnDists[i] = distsFromCurrVec[i];
                    nnIndices[i] = windowMaxSize - 1;

                    sumAddNewdists = sumAddNewdists + nnDists[i];
                }
            }

        }

        // Case 3: The partition the last point was deleted from is the same as the target partition of the new point to be added to,
        // and the i-th point belongs to that partition (i.e. at least one point remained in the partition after the deletion and before
        // the addition).
        else if (labelToBeDeleted == targetPartition && partitionLabels[i] == targetPartition)
        {

            // If only one point remained in the partition after the deletion and before the addition:
            if ( numPointsPartn[targetPartition] == 1 )
            {

                // Update the NN statistics for the "new" two-point partition.
                nnIndices[i] = windowMaxSize - 1;
                nnIndices.push_back(i);

                nnDists[i] = distsFromCurrVec[i];
                nnDists.push_back(nnDists[i]);

                avgNNDistPartitions[targetPartition] = nnDists[i];
                numPointsPartn[targetPartition] = 2;
                twoPointPartition = 1;
            }

            else if ( distsFromCurrVec[i] < nnDists[i] )
            {
                distsFromCurrVecTP.push_back( distsFromCurrVec[i] );
                tpIndices.push_back(i);

                sumOldNNdists = sumOldNNdists + nnDists[i];

                nnDists[i] = distsFromCurrVec[i];
                nnIndices[i] = windowMaxSize - 1;

                sumNewNNdists = sumNewNNdists + nnDists[i];
            }

            else if ( nnIndices[i] == indexToBeDeleted )
            {
                distsFromCurrVecTP.push_back( distsFromCurrVec[i] );
                tpIndices.push_back(i);

                sumOldNNdists = sumOldNNdists + nnDists[i];

                std::vector<double> memberDistsFromVect;
                std::vector<int> memberIndices;

                for(unsigned int j = 0; j < windowMaxSize; j++)
                {
                    if ( partitionLabels[j] == labelToBeDeleted )
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
                nnIndices[i] = memberIndices[tempIdx];

                auto newNNdistFromVect = *std::min_element( memberDistsFromVect.begin(), memberDistsFromVect.end() );
                nnDists[i] = newNNdistFromVect;
                sumNewNNdists = sumNewNNdists + newNNdistFromVect;
            }

            else
            {
                distsFromCurrVecTP.push_back( distsFromCurrVec[i] );
                tpIndices.push_back(i);
            }
        }

    }

    std::vector<double> distMatLastRow(windowMaxSize);  // The last row of the upper triangular distance matrix is a vector of 0s.
    // std::vector<double> distMatLastRow = distsFromCurrVec;
    // distMatLastRow.push_back(0);
    //pPack->complex->distMatrix->push_back( distMatLastRow );
    distMatrix.push_back( distMatLastRow );

    pPack.complex->setDistanceMatrix(&distMatrix);


    // Update the average NN distance of the partition from which the last point was deleted and of the one to which the new point
    // is being added.
    // If the new point is being added to a partition that does not exist in the window (after the last deletion):
    if ( avgNNDistPartitions.count( targetPartition ) == 0 )
    {
        avgNNDistPartitions[targetPartition] = -1;
        nnIndices.push_back(-1);
        nnDists.push_back(-1);
        numPointsPartn[targetPartition] = 1;

    }
    if ( labelToBeDeleted != targetPartition )
    {

        // Update the average NN distance of the partition from which the last point was deleted.
        if ( numPointsPartn[labelToBeDeleted] > 1 )
        {

            int n = numPointsPartn[labelToBeDeleted];
            double avgDel = avgNNDistPartitions[labelToBeDeleted];
            avgNNDistPartitions[labelToBeDeleted] = ( (n+1)*avgDel - nnDistToBeDeleted - sumDelOldNNdists + sumDelNewdists ) / n;
        }

        // Update the average NN distance of the partition to which the new point is being added.
        if ( numPointsPartn[targetPartition] >= 2 && twoPointPartition == 0 )
        {
            auto tempNNidx = std::min_element(distsFromCurrVecTP.begin(), distsFromCurrVecTP.end()) - distsFromCurrVecTP.begin();
            nnIndices.push_back( tpIndices[tempNNidx] );

            auto nnDistFromNewVect = *std::min_element( distsFromCurrVecTP.begin(), distsFromCurrVecTP.end() );
            nnDists.push_back( nnDistFromNewVect );

            int n = numPointsPartn[targetPartition];
            double avgAdd = avgNNDistPartitions[targetPartition];
            avgNNDistPartitions[targetPartition] = ( n*avgAdd + nnDistFromNewVect - sumAddOldNNdists + sumAddNewdists ) / (n+1);

            numPointsPartn[targetPartition] = n + 1;

        }
    }
    else if ( labelToBeDeleted == targetPartition
              && numPointsPartn[targetPartition] >= 2
              && twoPointPartition == 0 )
    {
        auto tempNNidx = std::min_element(distsFromCurrVecTP.begin(), distsFromCurrVecTP.end()) - distsFromCurrVecTP.begin();
        nnIndices.push_back( tpIndices[tempNNidx] );

        auto nnDistFromNewVect = *std::min_element( distsFromCurrVecTP.begin(), distsFromCurrVecTP.end() );
        nnDists.push_back( nnDistFromNewVect );

        int n = numPointsPartn[targetPartition] + 1;
        double avgNNd = avgNNDistPartitions[targetPartition];

        avgNNDistPartitions[targetPartition] = ( n*avgNNd - nnDistToBeDeleted - sumOldNNdists + sumNewNNdists + nnDistFromNewVect ) / n;
        numPointsPartn[targetPartition] = n;
    }

    return;
}

template<typename nodeType>
bool slidingWindow<nodeType>::nnBasedEvaluator(std::vector<double>& currentVector, std::vector<std::vector<double>>& windowValues)
{
    utils ut;
    float f1{ 4 };
    float f2{ 0.25 };
    float f3{ 5 };
    distsFromCurrVec.clear();

    // Compute the distances from the current vector to the existing ones in the window.
    distsFromCurrVec = ut.nearestNeighbors(currentVector, windowValues);

    // Find the distance from the current vector to its nearest neighbor in the window.
    auto nnDistCurrVec = *std::min_element( distsFromCurrVec.begin(), distsFromCurrVec.end() );

    // std::cout << nnDistCurrVec << std::endl;

    if (avgNNDistPartitions.size() == 1)  // If the window is 'pure':
    {

        if (nnDistCurrVec == 0) {

//            std::cout << "1. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"  << '\n';
            std::cout << "pointCounter = " << pointCounter << '\n';
//            for (auto const& pair: avgNNDistPartitions) {
//                std::cout << "{" << pair.first << ": " << pair.second << "}\n";
//            }
            std::cout << "1 ==============================================================================================" << '\n';

            return false;
        }

        // Find the average nearest neighbor distance in the single 'partition' in the window.
        int currentLabel = partitionLabels[0];
        // std::cout << "currentLabel: " << currentLabel << '\n';

        auto avgNNDistSinglePartition = avgNNDistPartitions[currentLabel];

        // std::cout << "avgNNDistSinglePartition: " << avgNNDistSinglePartition << std::endl;

        if (avgNNDistSinglePartition <= f2 && nnDistCurrVec <= 1) {
//            std::cout << "2. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"  << '\n';
//
            std::cout << "pointCounter = " << pointCounter << '\n';
//            for (auto const& pair: avgNNDistPartitions) {
//                std::cout << "{" << pair.first << ": " << pair.second << "}\n";
//            }
            std::cout << "2 ==============================================================================================" << '\n';

            return false;
        }


        if (avgNNDistSinglePartition == 0 || nnDistCurrVec / avgNNDistSinglePartition > f1)
        {
            // std::cout << "3. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"  << '\n';
            // In this case, the oldest point will be deleted from the window.
            labelToBeDeleted = partitionLabels[0];
            indexToBeDeleted = 0;

            // Since the window was 'pure', the partition label of the new point to be added is (existing label + 1).
            targetPartition = labelToBeDeleted + 1;

            deleteNNstats();

            // Delete the vector from the front of the sliding window.
            windowValues.erase( windowValues.begin() );

            // Update the distance matrix and the NN statistics.
            updateStats();

            std::cout << "1. Point Added" << '\n';

            std::cout << "pointCounter = " << pointCounter << '\n';
//            for (auto const& pair: avgNNDistPartitions) {
//                std::cout << "{" << pair.first << ": " << pair.second << "}\n";
//            }
            std::cout << "==============================================================================================" << '\n';

            return true;
        }

    }

    else {   // If the window is NOT "pure":

        // Determine the partition membership of the current vector. In particular, check if the current vector can be assigned to
        // its nearest partition. If it cannot be assigned to its nearest partition, create a new partition with only the current vector.

        // Find the partition that is nearest to the current vector.
        auto nearestPartnIdx = std::min_element(distsFromCurrVec.begin(), distsFromCurrVec.end()) - distsFromCurrVec.begin();
        int nearestPartition = partitionLabels[nearestPartnIdx];

        // Find the max. of the existing partition labels.
        int maxLabel = *std::max_element( partitionLabels.begin(), partitionLabels.end() );

        // If the (nearest neighbor) distance from the current vector to any partition is 0, assign the current vector to that partition.
        if ( nnDistCurrVec == 0 ) {
           targetPartition = nearestPartition;
        }

        // If the minimum distance from the current vector to existing partitions is higher than f3, assign a new partition to the current vector.
        else if ( nnDistCurrVec > f3 ) {
            targetPartition = maxLabel + 1;
        }
        else if ( avgNNDistPartitions[nearestPartition] == -1 && nnDistCurrVec <= f3 ) {
            targetPartition = nearestPartition;
        }
        else if ( avgNNDistPartitions[nearestPartition] <= f2 && nnDistCurrVec <= 1 ) {
            targetPartition = nearestPartition;
        }
        else if ( avgNNDistPartitions[nearestPartition] > 0 && nnDistCurrVec / avgNNDistPartitions[nearestPartition] <= f1 ) {
            targetPartition = nearestPartition;
        }
        else {
            targetPartition = maxLabel + 1;
        }

        // Determine which point would be deleted from the sliding window.
        // A partition is considered outdated if it did not receive any new point for more than the last 25 insertions.
        int timeToBeOutdated{ 25 };

        // First check if there are outdated partitions. If there are, find the smallest one of them.
        int smallestOutdated{ -1 };   // This will store the label of the smallest outdated partition.
        int numPtsSmallestOutdated{ -1 };  // This will store the number of points in the smallest outdated partition.

        for(auto op : maxKeys) {
            if ( (key - op.second) > timeToBeOutdated ) {   // If op is outdated:
                if ( numPtsSmallestOutdated == -1 ) {   // If this is the first time an outdated partition is encountered:
                    numPtsSmallestOutdated = numPointsPartn[op.first];
                    smallestOutdated = op.first;
                }
                else if ( numPointsPartn[op.first] < numPtsSmallestOutdated ) {   // If a smaller outdated partition is encountered:
                    // Update our records for the smallest outdated partition.
                    numPtsSmallestOutdated = numPointsPartn[op.first];
                    smallestOutdated = op.first;
                }
            }
        }

        if (smallestOutdated != -1) {   // If there is (are) outdated partition(s):

            // std::cout << "4. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"  << '\n';
            // Delete the oldest point of the smallest outdated partition (and its associated statistics) from the sliding window.
            labelToBeDeleted = smallestOutdated;
            indexToBeDeleted = std::find( partitionLabels.begin(), partitionLabels.end(), smallestOutdated ) - partitionLabels.begin();

            deleteNNstats();

            windowValues.erase( windowValues.begin() + indexToBeDeleted );

            // Update the distance matrix and the NN statistics.
            updateStats();

            std::cout << "2. Point Added" << '\n';

            std::cout << "pointCounter = " << pointCounter << '\n';
//            for (auto const& pair: avgNNDistPartitions) {
//                std::cout << "{" << pair.first << ": " << pair.second << "}\n";
//            }
            std::cout << "==============================================================================================" << '\n';

        }

        else {   // There is no outdated partition in the window:
            if ( targetPartition != nearestPartition ) {   // If the current vector was assigned a new partition:

                // std::cout << "5. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"  << '\n';
                // In this case, the oldest point will be deleted from the window.
                labelToBeDeleted = partitionLabels[0];
                indexToBeDeleted = 0;

                deleteNNstats();

                // Delete the vector from the front of the sliding window.
                windowValues.erase( windowValues.begin() );

                // Update the distance matrix and the NN statistics.
                updateStats();

                std::cout << "3. Point Added" << '\n';

                std::cout << "pointCounter = " << pointCounter << '\n';
//                for (auto const& pair: avgNNDistPartitions)
//                {
//                    std::cout << "{" << pair.first << ": " << pair.second << "}\n";
//                }
                std::cout << "==============================================================================================" << '\n';

            }
            else {   // The current vector is assigned to one of the existing partitions:
                // In this case, make sure the point to be deleted does not belong to the target partition. In particular,
                // we'll delete the oldest point from a partition label != targetPartition label.

                // std::cout << "6. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"  << '\n';

                // Find the first occurrence of a partition label != targetPartition label.
                for (auto delIndex = 0; delIndex < partitionLabels.size(); delIndex++) {
                    if ( partitionLabels[delIndex] != targetPartition ) {
                        indexToBeDeleted = delIndex;
                        break;
                    }

                }

                // std::cout << "indexToBeDeleted = " << indexToBeDeleted << '\n';
                // indexToBeDeleted = std::find( partitionLabels.begin(), partitionLabels.end(), !targetPartition ) - partitionLabels.begin();
                labelToBeDeleted = partitionLabels[indexToBeDeleted];

                deleteNNstats();
                windowValues.erase( windowValues.begin() + indexToBeDeleted );

                // Update the distance matrix and the NN statistics.
                updateStats();

                std::cout << "4. Point Added" << '\n';

                std::cout << "pointCounter = " << pointCounter << '\n';
//                for (auto const& pair: avgNNDistPartitions)
//                {
//                    std::cout << "{" << pair.first << ": " << pair.second << "}\n";
//                }
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
template<typename nodeType>
pipePacket<nodeType> slidingWindow<nodeType>::runPipe(pipePacket<nodeType> inData)
{
    readInput rp;
	pPack = inData;
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

        while(rp.streamRead(currentVector))
        {
            // Initialize the sliding window. During the initialization, let's assume all points from the stream
            // belong to Partition 0.
            if (windowValues.size() < windowMaxSize)
            {
                windowValues.push_back(currentVector);
                windowKeys.push_back(key);
                partitionLabels.push_back(targetPartition);
                key++;

                //If we've reached window max size, generate the initial complex
                if (windowValues.size() == windowMaxSize)
                {
                    std::cout << "Initializing complex" << std::endl;

                    inData.workData = windowValues;

                    distMatrix.resize(inData.workData.size(), std::vector<double>(inData.workData.size(),0));

					//Iterate through each vector
                    for(unsigned i = 0; i < inData.workData.size(); i++)
                    {
                        if(!inData.workData[i].empty())
                        {
                            std::vector<double> distsFromCurrVect;

                            for(unsigned j = 0; j < inData.workData.size(); j++)
                            {
                                if (j < i)
                                {
                                    distsFromCurrVect.push_back( distMatrix[j][i] );
                                }
                                else if (j > i)
                                {
                                    //Calculate vector distance
                                    auto dist = this->ut.vectors_distance(inData.workData[i], inData.workData[j]);
                                    distMatrix[i][j] = dist;
                                    distsFromCurrVect.push_back( dist );
                                }

                            }

                            auto tempIndex = std::min_element(distsFromCurrVect.begin(), distsFromCurrVect.end()) - distsFromCurrVect.begin();
                            if (tempIndex < i)
                                nnIndices.push_back(tempIndex);
                            else
                                nnIndices.push_back(tempIndex + 1);

                            auto nnDistFromCurrVect = *std::min_element( distsFromCurrVect.begin(), distsFromCurrVect.end() );
                            nnDists.push_back( nnDistFromCurrVect );
                        }
                    }

                    inData.complex->setDistanceMatrix(&distMatrix);

					for(auto a : windowValues)
						inData.complex->insert();

                    // Set the stream evaluator
                    inData.complex->setStreamEvaluator(&nnBasedEvaluator);

                    std::cout << "Returning from complex initializer" << std::endl;

                    // Find the average nearest neighbor distance in the existing partition (i.e. Partition 0).
                    auto avgNNDistPartition0 = std::accumulate(nnDists.begin(), nnDists.end(), 0.0) / nnDists.size();
                    avgNNDistPartitions[targetPartition] = avgNNDistPartition0;
                    numPointsPartn[targetPartition] = windowMaxSize;
                    maxKeys[targetPartition] = key - 1;
                }


            }
            else
            {
                if(inData.complex->insertIterative(currentVector, windowValues, keyToBeDeleted, indexToBeDeleted))
                {
                    // inData.complex->deleteIndexRecurse( keyToBeDeleted );

                    // Insert the current vector, its key and partition label into the rear ends of the corresponding containers.
                    windowValues.push_back(currentVector);
                    windowKeys.push_back(key);
                    partitionLabels.push_back(targetPartition);
                    maxKeys[targetPartition] = key;  // Insert or update the maxKey of the target partition.
                    key++;

                }
            }

            //Check if we've gone through 100 points
            if(pointCounter % 500 == 0 && pointCounter >= windowMaxSize)
            {
                // Build and trigger remaining pipeline. It should only require the computation of persistence
                // intervals from the complex being maintained.
                std::cout << "pointCounter: " << pointCounter << "\tSimplex Count: " << inData.complex->simplexCount() << "\tVertex Count: " << inData.complex->vertexCount() << std::endl;
                std::cout << "\tWindowSize: " << windowValues.size() << std::endl;
                inData.workData = windowValues;
                runSubPipeline(inData);

            }
            pointCounter++;

            currentVector.clear();

        }
        //Probably want to trigger the remaining pipeline one last time...
        if((pointCounter - 1) % 500 != 0)
        {
            std::cout << "pointCounter: " << pointCounter << "\tSimplex Count: " << inData.complex->simplexCount() << "\tVertex Count: " << inData.complex->vertexCount() << std::endl;

            runSubPipeline(inData);
        }
        this->ut.writeLog("slidingWindow", "\tSuccessfully evaluated " + std::to_string(pointCounter) + " points");

        writeComplexStats(inData);
    }

    return inData;
}

template<typename nodeType>
void slidingWindow<nodeType>::writeComplexStats(pipePacket<nodeType> &inData)
{
    if(inData.complex->stats.size() > 30)
    {
        std::ofstream file ("output/complexStats.csv");

        file << inData.complex->stats << std::endl;

        file.close();

    }
}

template<typename nodeType>
void slidingWindow<nodeType>::runSubPipeline(pipePacket<nodeType> wrData)
{
    if(wrData.workData.size() == 0)
        return;

    pipePacket<nodeType> inData = wrData;
    outputData(inData);

	std::string pipeFuncts = "rips.fast";
	// std::string pipeFuncts = "rips";
    auto lim = count(pipeFuncts.begin(), pipeFuncts.end(), '.') + 1;
    subConfigMap["fn"] = "_" + std::to_string(repCounter);
    repCounter++;

    //For each '.' separated pipeline function (count of '.' + 1 -> lim)
    for(unsigned i = 0; i < lim; i++)
    {
        auto curFunct = pipeFuncts.substr(0,pipeFuncts.find('.'));
        pipeFuncts = pipeFuncts.substr(pipeFuncts.find('.') + 1);

        //Build the pipe component, configure and run
        auto cp = basePipe<nodeType>::newPipe(curFunct, "simplexTree");

        //Check if the pipe was created and configure
        if(cp != 0 && cp->configPipe(subConfigMap))
        {
            //Run the pipe function (wrapper)
            cp->runPipeWrapper(inData);
        }
        else
        {
            std::cout << "LHF : Failed to configure pipeline: " << curFunct << std::endl;
        }
    }


    return;
}

// configPipe -> configure the function settings of this pipeline segment
template<typename nodeType>
bool slidingWindow<nodeType>::configPipe(std::map<std::string, std::string> configMap)
{
    std::string strDebug;
    this->subConfigMap = configMap;

    auto pipe = configMap.find("debug");
    if(pipe != configMap.end())
    {
        this->debug = std::atoi(configMap["debug"].c_str());
        strDebug = configMap["debug"];
    }
    pipe = configMap.find("outputFile");
    if(pipe != configMap.end())
        this->outputFile = configMap["outputFile"].c_str();

    this->ut = utils(strDebug, this->outputFile);

    pipe = configMap.find("inputFile");
    if(pipe != configMap.end())
        this->inputFile = configMap["inputFile"].c_str();

    pipe = configMap.find("epsilon");
    if(pipe != configMap.end())
        this->epsilon = std::atof(configMap["epsilon"].c_str());
    else return false;

    pipe = configMap.find("dimensions");
    if(pipe != configMap.end())
    {
        this->dim = std::atoi(configMap["dimensions"].c_str());
    }
    else return false;

    this->configured = true;
    this->ut.writeDebug("slidingWindow","Configured with parameters { input: " + configMap["inputFile"] + ", dim: " + configMap["dimensions"] + ", eps: " + configMap["epsilon"] + ", debug: " + strDebug + ", outputFile: " + this->outputFile + " }");

    return true;
}

// outputData -> used for tracking each stage of the pipeline's data output without runtime
template<typename nodeType>
void slidingWindow<nodeType>::outputData(pipePacket<nodeType> inData)
{
    std::ofstream file ("output/" + this->pipeType + "_" + std::to_string(this->repCounter) + "_output.csv");

    for(auto a : inData.workData)
    {
        for(auto d : a)
        {
            file << d << ",";
        }
        file << "\n";
    }
    file << std::endl;

    file.close();
    return;
}

template class slidingWindow<simplexNode>;
template class slidingWindow<alphaNode>;
