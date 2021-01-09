/*
 *Interpolate and merge Betties for all diffrent executions
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
interpolateMergePipe::interpolateMergePipe(){
	pipeType = "interpolateMerge";
	return;
}


void distBuildSCExecutePHPipe::runPipe(pipePacket &inData){

std::vector<bettiBoundaryTableEntry> finalbettiTable;
std::ofstream file("bettisequence.csv");
for(auto entry : bTbs){
 //   file<<entry.bettiTab.size()<<","<<entry.numpts<<","<<entry.maxEpsilon;
 //   file<<endl;
    for(auto a : entry.bettiTab){
        file<<a.bettiDim<<","<<a.birth<<","<<a.death<<endl;
        bool found = false;
        for(auto x : finalbettiTable)
        {
           if(a.bettiDim==x.bettiDim&&a.birth==x.birth&&a.death==x.death&&a.boundaryPoints==x.boundaryPoints)
               {
                found = true;
                break;
               }
        }
        if(!found)
            finalbettiTable.push_back(a);
    }
 //   file<<endl;
}
file<<endl<<"Final Betti Table"<<endl;
for(auto a : finalbettiTable){
        file<<a.bettiDim<<","<<a.birth<<","<<a.death;
        for(auto k : a.boundaryPoints)
            file<<k<<",";
    file<<endl;
}
file.close();
}