 /*
 *Build Complexes for all the valid_d_simplexmeshes
 */

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <functional>
#include <set>
#include <algorithm>
#include "distMatrixPipe.hpp"
#include "utils.hpp"

// basePipe constructor
distBuildSCExecutePHPipe::distBuildSCExecutePHPipe(){
	pipeType = "distBuildSCExecutePH";
	return;
}

// Generic template
template <typename T>

pair<T, bool> distBuildSCExecutePHPipe:: getNthElement(set<T>& searchSet,unsigned index)
{
    pair<T, bool> result;

    // Check if index is valid or not
    if (searchSet.size() > index) {
        result.first
            = *(std::next(
                searchSet.begin(),
                index));
        result.second = true;
    }

    else
        result.second = false;

    // Return the pair
    return result;
}
	//Will embed this function into complex pipeline;

std::vector<set<simplexNode_P, cmpByWeight>> distBuildSCExecutePHPipe::buildValidSimplicialComplex(vector<set<unsigned>> dsimplexes,pipePacket &inData){

   vector<set<simplexNode_P, cmpByWeight>> dim;
   for(int i =0;i<=maxDimension;i++)
       dim.push_back({});


   for(auto simplex : dsimplexes)
   {
    unsigned int pow_set_size = pow(2, simplex.size());
    unsigned counter, j;
     for(counter = 1; counter < pow_set_size; counter++)
    {
        double weight =0;
        set<unsigned> gensimp;
    for(j = 0; j < simplex.size(); j++)
    {
       if(counter & (1 << j))
            {pair<unsigned, bool> index = getNthElement(simplex,j);
            if(index.second){
               unsigned indnew = index.first;
              for(auto x : gensimp)
                if(weight<inData.distMatrix[x][indnew])
                    weight = inData.distMatrix[x][indnew];
              gensimp.insert(indnew);
            }
            }

    }

     simplexNode_P tot = std::make_shared<simplexNode>(simplexNode(gensimp, weight));
     tot->hash1 = simplexHash(gensimp);
     dim[gensimp.size()-1].insert(tot);
    }
}

return dim;
}

std::vector<bettiBoundaryTableEntry> distBuildSCExecutePHPipe::computePersistence(){
	//Will embed this function into persistence pipeline;
}

// runPipe -> Run the configured functions of this pipeline segment
void distBuildSCExecutePHPipe::runPipe(pipePacket &inData){
	
for(auto validSC : inData.VUdsimplexes)
{
	std::vector<set<simplexNode_P, cmpByWeight>> simplexList;          \\ Complex 
	simplexList = buildValidSimplicialComplex(validSC,inData);
	std::vector<bettiBoundaryTableEntry> bettiTable = computePersistence(simplexList);
	bettiTables BT =  bettiTables(bettiTable,(*(simplexList[1].rbegin()))->weight,simplexList[0].size()) ;
	inData.bTbs.insert(BT);
}	
	return;
	
}

