 /*
 * enumerateUVSC hpp + cpp extend the basePipe class for enumerating the valid unique d_simplexes
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
enumerateUVSCPipe::enumerateUVSCPipe(){
	pipeType = "enumerateUVSC";
	return;
}
// runPipe -> Run the configured functions of this pipeline segment
void enumerateUVSCPipe::runPipe(pipePacket &inData){
std::vector<unsigned> permut(inData.inputData.size());
unsigned i = 0;            // incrementor
for_each(permut.begin(), permut.end(), [&](unsigned& item) { ++i; item += i;});
unsigned n = inData.inputData.size();
 do {
vector<set<unsigned>> intersectionfreemesh;
bool firstflag = true;


//Curently it is 2-simplexes need to generalize to d dimensions
for(unsigned i=0;i<inData.inputData.size();i++)
    for(unsigned j=i+1;j<inData.inputData.data.size();j++)
      for(unsigned k = j+1;k<indata.inputData.size();k++)
        {
          set<unsigned> newsimplex;
          newsimplex.insert(permut[i]);
          newsimplex.insert(permut[j]);
          newsimplex.insert(permut[k]);
          if(firstflag)
                intersectionfreemesh.push_back(newsimplex);
          vector<vector<double>> p;
          p.push_back(data[permut[i]]);
          p.push_back(data[permut[j]]);
          p.push_back(data[permut[k]]);

          bool intersection =false;
          if(!firstflag){

          unsigned pl =0;
                for(auto oldsimplex : intersectionfreemesh){

                    vector<vector<double>> q;
                    for(auto ind:oldsimplex)
                        q.push_back(data[ind]);

                     TriPoint t1[] = {TriPoint(p[0][0],p[0][1]),TriPoint(p[1][0],p[1][1]),TriPoint(p[2][0],p[2][1])};
	                 TriPoint t2[] = {TriPoint(q[0][0],q[0][1]),TriPoint(q[1][0],q[1][1]),TriPoint(q[2][0],q[2][1])};

                    if(d_simplexIntersection::TriTri2D(t1, t2)){
                        intersection = true;
                        break;
                    }
                    pl++;
                }

                if(!intersection)
                    intersectionfreemesh.push_back(newsimplex);
          }
          firstflag = false;
        }
bool present = false;
for(auto x : inData.VUdsimplexes){
    if(x==intersectionfreelist)
    {present = true;
     break;
    }
}
if(!present)
    inData.VUdsimplexes.push_back(intersectionfreemesh);

  } while (next_permutation(permut, permut + n));

	
	return;
}

