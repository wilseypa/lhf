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
#include "distBuildSCExecutePH.hpp"
#include "utils.hpp"

// basePipe constructor
distBuildSCExecutePHPipe::distBuildSCExecutePHPipe(){
	pipeType = "distBuildSCExecutePH";
	return;
}

std::vector<bettiBoundaryTableEntry>  vscPersistence::lightPersistence(){

	std::vector<bettiBoundaryTableEntry> bettiTable;

	return bettiTable;
}	// Complex 
// runPipe -> Run the configured functions of this pipeline segment
void distBuildSCExecutePHPipe::runPipe(pipePacket &inData){
	
for(auto validSC : inData.VUdsimplexes)
{
	std::vector<std::set<simplexNode_P, cmpByWeight>> simplicialComplex;          // Complex 
	simplicialComplex =  inData.complex->buildValidSimplicialComplex(validSC);
	
	double aEW=0;
        for(auto edge : simplicialComplex[1])
		aEW = aEW + (*edge).weight;
	aEW/= (double)simplicialComplex[1].size();

	auto vPH = vscPersistence(simplicialComplex);
        	
	std::vector<bettiBoundaryTableEntry> bettiTable = vPH.lightPersistence();
	

	bettiTables BT =  bettiTables(bettiTable,(*(simplicialComplex[1].rbegin()))->weight,aEW) ;

        inData.bTbs.insert(BT);
}	
	return;
	
}

