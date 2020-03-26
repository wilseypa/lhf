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
#include <functional>
#include "slidingWindow.hpp"
#include "utils.hpp"
#include "readInput.hpp"


// basePipe constructor
slidingWindow::slidingWindow()
{
    pipeType = "SlidingWindow";
    return;
}


bool nnBasedEvaluator(std::vector<double>& currentVector, std::vector<std::vector<double>>& windowValues, EvalParams& defaultVals)
{
    utils ut;
    float f1{ 4 };
    float f2{ 0.25 };
    float f3{ 5 };
    defaultVals.distsFromCurrVec.clear();

    // Compute the distances from the current vector to the existing ones in the window.
    defaultVals.distsFromCurrVec = ut.nearestNeighbors(currentVector, windowValues);

    // Find the distance from the current vector to its nearest neighbor in the window.
    auto nnDistCurrVec = *std::min_element( defaultVals.distsFromCurrVec.begin(), defaultVals.distsFromCurrVec.end() );

    if (defaultVals.avgNNDistPartitions.size() == 1)  // If the window is 'pure':
    {
        if (nnDistCurrVec == 0)
            return false;

        // Find the average nearest neighbor distance in the single 'partition' in the window.
        int currentLabel = defaultVals.partitionLabels[0];
        auto avgNNDistSinglePartition = defaultVals.avgNNDistPartitions[currentLabel];

        if (avgNNDistSinglePartition <= f2 && nnDistCurrVec <= 1)
            return false;

        if (avgNNDistSinglePartition == 0 || nnDistCurrVec / avgNNDistSinglePartition > f1)
        {
            // Delete the key (the lowest key) from the front of the list.
            defaultVals.keyToBeDeleted = defaultVals.windowKeys[0];
            defaultVals.windowKeys.erase( defaultVals.windowKeys.begin() );

            // Delete the label from the front of the list.
            defaultVals.labelToBeDeleted = defaultVals.partitionLabels[0];
            defaultVals.partitionLabels.erase( defaultVals.partitionLabels.begin() );

            defaultVals.indexToBeDeleted = 0; // In this case, the oldest point will be deleted from the window.

            defaultVals.nnIndices.erase( defaultVals.nnIndices.begin() );

            defaultVals.nnDistToBeDeleted = defaultVals.nnDists[0];
		    defaultVals.nnDists.erase( defaultVals.nnDists.begin() );

            // Delete the corresponding distance value from the list of distances from the current vector
            // to the existing ones in the window.
            defaultVals.distsFromCurrVec.erase( defaultVals.distsFromCurrVec.begin() );

            // Delete the vector from the front of the sliding window.
            windowValues.erase( windowValues.begin() );

            // Since the window was 'pure', the partition label of the new point to be added is (existing label + 1).
            defaultVals.targetPartition = defaultVals.labelToBeDeleted + 1;

            // Decrement the number of points in the existing partition by 1.
            defaultVals.numPointsPartn[defaultVals.labelToBeDeleted] = defaultVals.windowMaxSize - 1;

            return true;
        }

    }

    else {   // If the window is NOT "pure":

//        // Create a dictionary to store the nearest neighbor distance from the current vector to each partition in the window.
//        std::unordered_map<int, double> nnDistsFrmCurrVecToPartns;
//        for(unsigned int i = 0; i < defaultVals.windowMaxSize; i++) {
//
//            // If the partition label of the i-th point does not already exist in nnDistsFrmCurrVecToPartns:
//            if ( nnDistsFrmCurrVecToPartns.count(defaultVals.partitionLabels[i]) == 0 ) {
//                nnDistsFrmCurrVecToPartns[defaultVals.partitionLabels[i]] = defaultVals.distsFromCurrVec[i];
//            }
//            else if ( defaultVals.distsFromCurrVec[i] < nnDistsFrmCurrVecToPartns[defaultVals.partitionLabels[i]] ) {
//                nnDistsFrmCurrVecToPartns[defaultVals.partitionLabels[i]] = defaultVals.distsFromCurrVec[i];
//            }
//        }



        // Determine the partition membership of the current vector. In particular, check if the current vector can be assigned to
        // its nearest partition. If it cannot be assigned to its nearest partition, create a new partition with only the current vector.

        // Find the partition that is nearest to the current vector.
        auto nearestPartnIdx = std::min_element(defaultVals.distsFromCurrVec.begin(), defaultVals.distsFromCurrVec.end()) - defaultVals.distsFromCurrVec.begin();
        int nearestPartition = defaultVals.partitionLabels[nearestPartnIdx];

        // Find the max. of the existing partition labels.
        int maxLabel = *std::max_element( defaultVals.partitionLabels.begin(), defaultVals.partitionLabels.end() );

        // If the (nearest neighbor) distance from the current vector to any partition is 0, assign the current vector to that partition.
        if ( nnDistCurrVec == 0 ) {
           defaultVals.targetPartition = nearestPartition;
        }

        // If the minimum distance from the current vector to existing partitions is higher than f3, assign a new partition to the current vector.
        else if ( nnDistCurrVec > f3 ) {
            defaultVals.targetPartition = maxLabel + 1;
        }
        else if ( defaultVals.avgNNDistPartitions[nearestPartition] == -1 && nnDistCurrVec <= f3 ) {
            defaultVals.targetPartition = nearestPartition;
        }
        else if ( defaultVals.avgNNDistPartitions[nearestPartition] <= f2 && nnDistCurrVec <= 1 ) {
            defaultVals.targetPartition = nearestPartition;
        }
        else if ( defaultVals.avgNNDistPartitions[nearestPartition] > 0 && nnDistCurrVec / defaultVals.avgNNDistPartitions[nearestPartition] <= f1 ) {
            defaultVals.targetPartition = nearestPartition;
        }
        else {
            defaultVals.targetPartition = maxLabel + 1;
        }

        // Determine which point would be deleted from the sliding window.
        // A partition is considered outdated if it did not receive any new point for more than the last 25 insertions.
        int timeToBeOutdated{ 25 };

        // First check if there are outdated partitions. If there are, find the smallest one of them.
        int smallestOutdated{ -1 };   // This will store the label of the smallest outdated partition.
        int numPtsSmallestOutdated{ -1 };  // This will store the number of points in the smallest outdated partition.

        for(auto op : defaultVals.maxKeys) {
            if ( (defaultVals.key - op.second) > timeToBeOutdated ) {   // If op is outdated:
                if ( numPtsSmallestOutdated == -1 ) {   // If this is the first time an outdated partition is encountered:
                    numPtsSmallestOutdated = defaultVals.numPointsPartn[op.first];
                    smallestOutdated = op.first;
                }
                else if ( defaultVals.numPointsPartn[op.first] < numPtsSmallestOutdated ) {   // If a smaller outdated partition is encountered:
                    // Update our records for the smallest outdated partition.
                    numPtsSmallestOutdated = defaultVals.numPointsPartn[op.first];
                    smallestOutdated = op.first;
                }
            }
        }

        if (smallestOutdated != -1) {   // If there is (are) outdated partition(s):

            // Delete the oldest point of the smallest outdated partition (and its associated statistics) from the sliding window.
            defaultVals.labelToBeDeleted = smallestOutdated;
            defaultVals.indexToBeDeleted = std::find( defaultVals.partitionLabels.begin(), defaultVals.partitionLabels.end(), smallestOutdated ) - defaultVals.partitionLabels.begin();
            defaultVals.keyToBeDeleted = defaultVals.windowKeys[defaultVals.indexToBeDeleted];

            defaultVals.windowKeys.erase( defaultVals.windowKeys.begin() + defaultVals.indexToBeDeleted );
            defaultVals.partitionLabels.erase( defaultVals.partitionLabels.begin() + defaultVals.indexToBeDeleted );

            defaultVals.nnIndices.erase( defaultVals.nnIndices.begin() + defaultVals.indexToBeDeleted );

            defaultVals.nnDistToBeDeleted = defaultVals.nnDists[defaultVals.indexToBeDeleted];
            defaultVals.nnDists.erase( defaultVals.nnDists.begin() + defaultVals.indexToBeDeleted );

            defaultVals.distsFromCurrVec.erase( defaultVals.distsFromCurrVec.begin() + defaultVals.indexToBeDeleted );

            windowValues.erase( windowValues.begin() + defaultVals.indexToBeDeleted );

            defaultVals.numPointsPartn[defaultVals.labelToBeDeleted] = defaultVals.numPointsPartn[defaultVals.labelToBeDeleted] - 1;

            // If there are no more points left in the partition from which the deletion took place:
            if (defaultVals.numPointsPartn[defaultVals.labelToBeDeleted] == 0)
            {
                defaultVals.avgNNDistPartitions.erase(defaultVals.labelToBeDeleted);
                defaultVals.numPointsPartn.erase(defaultVals.labelToBeDeleted);
                defaultVals.maxKeys.erase(defaultVals.labelToBeDeleted);
            }

        }

        return true;
    }

    return false;
}

// runPipe -> Run the configured functions of this pipeline segment
pipePacket slidingWindow::runPipe(pipePacket inData)
{
    utils ut;
    readInput rp;

    // For this pipe, we construct a sub-pipeline:
    //		1. Read data vector by vector, push into slidingWindow evaluation
    //		2. IF the point is to be inserted, push into the simplexTree (LRU ordered)
    //		3. IF 100 points have passed, generate persistence intervals
    //
    // Loop this subpipeline until there's no more data

    EvalParams defaultVals{ 200, 0, 0 };

    std::vector<std::vector<double>> windowValues;

    std::vector<double> currentVector;
    if(rp.streamInit(inputFile))
    {
        int pointCounter = 1;

        while(rp.streamRead(currentVector))
        {
            // Initialize the sliding window. During the initialization, let's assume all points from the stream
            // belong to Partition 0.
            if (windowValues.size() < defaultVals.windowMaxSize)
            {
                windowValues.push_back(currentVector);
                defaultVals.windowKeys.push_back(defaultVals.key);
                defaultVals.partitionLabels.push_back(defaultVals.targetPartition);
                defaultVals.key++;

                //If we've reached window max size, generate the initial complex
                if (windowValues.size() == defaultVals.windowMaxSize)
                {
                    std::cout << "Initializing complex" << std::endl;

                    inData.originalData = windowValues;
                    runComplexInitializer(inData, defaultVals.nnIndices, defaultVals.nnDists);

                    // Set the stream evaluator
                    inData.complex->setStreamEvaluator(&nnBasedEvaluator);

                    std::cout << "Returning from complex initializer" << std::endl;

                    // Find the average nearest neighbor distance in the existing partition (i.e. Partition 0).
                    auto avgNNDistPartition0 = std::accumulate(defaultVals.nnDists.begin(), defaultVals.nnDists.end(), 0.0) / defaultVals.nnDists.size();
                    defaultVals.avgNNDistPartitions[defaultVals.targetPartition] = avgNNDistPartition0;
                    defaultVals.numPointsPartn[defaultVals.targetPartition] = defaultVals.windowMaxSize;
                    defaultVals.maxKeys[defaultVals.targetPartition] = defaultVals.key - 1;
                }

            }
            else
            {
                if(inData.complex->insertIterative(currentVector, windowValues, defaultVals))
                {

                    // Insert the current vector, its key and partition label into the rear ends of the corresponding containers.
                    windowValues.push_back(currentVector);
                    defaultVals.windowKeys.push_back(defaultVals.key);
                    defaultVals.partitionLabels.push_back(defaultVals.targetPartition);
                    defaultVals.maxKeys[defaultVals.targetPartition] = defaultVals.key;  // Insert or update the maxKey of the target partition.
                    defaultVals.key++;

                }
            }

            //Check if we've gone through 100 points
            if(pointCounter % 100 == 0 && pointCounter > windowMaxSize)
            {
                // Build and trigger remaining pipeline. It should only require the computation of persistence
                // intervals from the complex being maintained.
                std::cout << "pointCounter: " << pointCounter << "\tSimplex Count: " << inData.complex->simplexCount() << "\tVertex Count: " << inData.complex->vertexCount() << std::endl;
                std::cout << "\tWindowSize: " << windowValues.size() << std::endl;
                inData.originalData = windowValues;
                runSubPipeline(inData);

            }
            pointCounter++;

            currentVector.clear();

        }
        //Probably want to trigger the remaining pipeline one last time...
        if((pointCounter - 1) % 100 != 0)
        {
            std::cout << "pointCounter: " << pointCounter << "\tSimplex Count: " << inData.complex->simplexCount() << "\tVertex Count: " << inData.complex->vertexCount() << std::endl;

            runSubPipeline(inData);
        }
        ut.writeLog("slidingWindow", "\tSuccessfully evaluated " + std::to_string(pointCounter) + " points");

        writeComplexStats(inData);
    }

    return inData;
}

void slidingWindow::writeComplexStats(pipePacket &inData)
{
    if(inData.complex->stats.size() > 30)
    {
        std::ofstream file ("output/complexStats.csv");

        file << inData.complex->stats << std::endl;

        file.close();

    }
}

void slidingWindow::runSubPipeline(pipePacket wrData)
{
    if(wrData.originalData.size() == 0)
        return;

    pipePacket inData = wrData;
    outputData(inData);

    std::string pipeFuncts = "rips.persistence";
    auto lim = count(pipeFuncts.begin(), pipeFuncts.end(), '.') + 1;
    subConfigMap["fn"] = "_" + std::to_string(repCounter);
    repCounter++;

    //For each '.' separated pipeline function (count of '.' + 1 -> lim)
    for(unsigned i = 0; i < lim; i++)
    {
        auto curFunct = pipeFuncts.substr(0,pipeFuncts.find('.'));
        pipeFuncts = pipeFuncts.substr(pipeFuncts.find('.') + 1);

        //Build the pipe component, configure and run
        auto *bp = new basePipe();
        auto *cp = bp->newPipe(curFunct, "simplexTree");

        //Check if the pipe was created and configure
        if(cp != 0 && cp->configPipe(subConfigMap))
        {
            //Run the pipe function (wrapper)
            inData = cp->runPipeWrapper(inData);
        }
        else
        {
            std::cout << "LHF : Failed to configure pipeline: " << curFunct << std::endl;
        }
    }


    return;
}

void slidingWindow::runComplexInitializer(pipePacket &inData, std::vector<int> &nnIndices, std::vector<double> &nnDists)
{
    //Initialize the complex and build other structures for maintaining NN, etc.
    //
    //	We need to exit this function by covering the distMatrix and neighGraph
    //		pipe functions



    //	1.	Create distance matrix (and compute other info)-------------

    utils ut;

    //Store our distance matrix
    std::vector<std::vector<double>> distMatrix (inData.originalData.size(), std::vector<double>(inData.originalData.size(),0));

    //Iterate through each vector
    for(unsigned i = 0; i < inData.originalData.size(); i++)
    {
        if(!inData.originalData[i].empty())
        {
            std::vector<double> distsFromCurrVect;

            for(unsigned j = 0; j < inData.originalData.size(); j++)
            {
                if (j < i)
                {
                    distsFromCurrVect.push_back( distMatrix[j][i] );
                }
                else if (j > i)
                {
                    //Calculate vector distance
                    auto dist = ut.vectors_distance(inData.originalData[i], inData.originalData[j]);
                    if(dist < epsilon)
                        inData.weights.insert(dist);
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

    inData.complex->setDistanceMatrix(distMatrix);

    inData.weights.insert(0.0);
    inData.weights.insert(epsilon);
    //std::sort(inData.weights.begin(), inData.weights.end(), std::greater<>());

    //------------------------------------------------------------------



    // 2. Insert into complex (build neighborhood graph) ---------------

    //Iterate through each vector, inserting into simplex storage
    for(unsigned i = 0; i < inData.originalData.size(); i++)
    {
        if(!inData.originalData[i].empty())
        {
            //insert data into the complex (SimplexArrayList, SimplexTree)
            inData.complex->insert(inData.originalData[i]);
        }
    }

    //------------------------------------------------------------------


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
void slidingWindow::outputData(pipePacket inData)
{
    std::ofstream file ("output/" + pipeType + "_" + std::to_string(repCounter) + "_output.csv");

    for(auto a : inData.originalData)
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
