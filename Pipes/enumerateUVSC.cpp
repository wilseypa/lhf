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
#include "enumerateUVSC.hpp"
#include "utils.hpp"

// basePipe constructor
enumerateUVSCPipe::enumerateUVSCPipe(){
	pipeType = "enumerateUVSC";
	return;
}
// runPipe -> Run the configured functions of this pipeline segment
void enumerateUVSCPipe::runPipe(pipePacket &inData){

// Doing permutaion for large point cloud is not feasible !! Hence optimization or approximation is required
std::vector<unsigned> permut(inData.inputData.size());
unsigned convex_hull_points = 14;
unsigned totaltriangles = 2*inData.inputData.size()-convex_hull_points-2;

/*
unsigned i = -1;            // incrementor
for_each(permut.begin(), permut.end(), [&](unsigned& item) { ++i; item += i;});
unsigned counter =0;
unsigned n = inData.inputData.size();
*/
std::vector<std::set<unsigned>> intersectionfreemesh;
std::vector<unsigned> intersectioncount(totaltriangles,0);

std::set<unsigned> newsimplex;
bool firstflag=true;;

std::vector<std::vector<double>> p;
std::vector<std::vector<double>> q;

int pk;
/* do {
	intersectionfreemesh.clear();
	firstflag = true;

//Curently it is 2-simplexes need to be generalized to d dimensions
*/
std::set<simplexNode_P,cmpByavgwt> dsimplexes;
std::vector<std::vector<double>> *distM = inData.complex->distMatrix; 
unsigned counter = 0,collin=0,noncollin=0;
for(unsigned i=0;i<inData.inputData.size();i++)
    for(unsigned j=i+1;j<inData.inputData.size();j++)
      for(unsigned k = j+1;k<inData.inputData.size();k++)
        { 
	  newsimplex.clear();
	  double a = (*distM)[i][j];
	  double b = (*distM)[i][k];
	  double c = (*distM)[j][k];
	  double s = (a+b+c)/2;
	  newsimplex.insert(i);
	  newsimplex.insert(j);
	  newsimplex.insert(k);
	  double maxWeight = std::max(a,std::max(b,c));
	  double perimeter = a+b+c;

	  simplexNode_P tot = std::make_shared<simplexNode>(newsimplex, maxWeight,perimeter);
  	  dsimplexes.insert(tot);
	  double area = sqrt(s*(s-a)*(s-b)*(s-c));
	  if(area==0)
               collin++;
	  else{
	       noncollin++;
               dsimplexes.insert(tot);
	  }
	  counter++;
	 
if(counter%1000000==0){
std::cout<<counter<<" "<<collin<<" "<<noncollin<<std::endl;
std::cout<<dsimplexes.size()<<" "<<(*dsimplexes.begin())->avgwt<<" "<<(*dsimplexes.rbegin())->avgwt<<std::endl;
}


	}

std::cout<<counter<<" "<<collin<<" "<<noncollin<<std::endl;
std::cout<<dsimplexes.size()<<" "<<(*dsimplexes.begin())->avgwt<<" "<<(*dsimplexes.rbegin())->avgwt<<std::endl;
counter = 0;
    for(auto x : dsimplexes){
	    
	  newsimplex.clear();
	  for(auto y : x->simplex)
         	 newsimplex.insert(y);
          if(firstflag)
                intersectionfreemesh.push_back(newsimplex);
	  p.clear();
          for(auto ind : newsimplex)
                p.push_back(inData.inputData[ind]);
          bool intersection =false;
          if(!firstflag){
		  unsigned intcounter = 0;
                for(auto oldsimplex : intersectionfreemesh){
   		    q.clear();
                    for(auto ind:oldsimplex)
                        q.push_back(inData.inputData[ind]);

                     TriPoint t1[] = {TriPoint(p[0][0],p[0][1]),TriPoint(p[1][0],p[1][1]),TriPoint(p[2][0],p[2][1])};
	             TriPoint t2[] = {TriPoint(q[0][0],q[0][1]),TriPoint(q[1][0],q[1][1]),TriPoint(q[2][0],q[2][1])};

                    if(d_simplexIntersection::TriTri2D(t1, t2)){
                        intersection = true;
			break;
                        intersectioncount[intcounter]++;
                       
                    }
		    intcounter++;
                }

                if(!intersection){
                    intersectionfreemesh.push_back(newsimplex);
		}
         }
	  if(totaltriangles == intersectionfreemesh.size())
	      break;
	  if(counter%100000==0)
	  {std::cout<<counter<<"  "<<totaltriangles<<" "<<intersectionfreemesh.size()<<std::endl;
	/*	 unsigned maxintersect = 0,index=0,cnt=0; 
   	 	 for(auto x:intersectioncount){
	   //      	std::cout<<" "<<x<<"-->";

		  if(x>maxintersect){
		    maxintersect = x;
		  index=cnt;}
		    cnt++;
		 }
	 std::cout<<maxintersect<<" "<<index<<" ";
	 intersectionfreemesh.erase(intersectionfreemesh.begin()+index);
	 std::fill(intersectioncount.begin(),intersectioncount.end(),0);
	 */
	  /*  int tp;
	    std::cin>>tp;
	    if(tp)
		    break;
		    */
		 
	  }
	  firstflag = false;
      counter ++;
    }
    inData.VUdsimplexes.insert(intersectionfreemesh);

/*	  
	  newsimplex.clear();
          newsimplex.insert(permut[i]);
          newsimplex.insert(permut[j]);
          newsimplex.insert(permut[k]);
          if(firstflag)
                intersectionfreemesh.push_back(newsimplex);
	  p.clear();
          p.push_back(inData.inputData[permut[i]]);
          p.push_back(inData.inputData[permut[j]]);
          p.push_back(inData.inputData[permut[k]]);

          bool intersection =false;
          if(!firstflag){

          unsigned pl =0;
                for(auto oldsimplex : intersectionfreemesh){
                    q.clear();
                    for(auto ind:oldsimplex)
                        q.push_back(inData.inputData[ind]);

                     TriPoint t1[] = {TriPoint(p[0][0],p[0][1]),TriPoint(p[1][0],p[1][1]),TriPoint(p[2][0],p[2][1])};
	             TriPoint t2[] = {TriPoint(q[0][0],q[0][1]),TriPoint(q[1][0],q[1][1]),TriPoint(q[2][0],q[2][1])};

                    if(d_simplexIntersection::TriTri2D(t1, t2)){
                        intersection = true;
                        break;
                    }
                    pl++;
                }

                if(!intersection){
                    intersectionfreemesh.push_back(newsimplex);
       		  }
	  }
          firstflag = false;


        }
    inData.VUdsimplexes.insert(intersectionfreemesh);
    counter++;
    if(counter%10000 ==0){	  
	    std::cout<<inData.VUdsimplexes.size()<<"-->";
	   // std::cin>>pk;
    }


  } while (next_permutation(permut.begin(), permut.end()));
*/
	return;
}


// configPipe -> configure the function settings of this pipeline segment
bool enumerateUVSCPipe::configPipe(std::map<std::string, std::string> &configMap){
	std::string strDebug;
	
	auto pipe = configMap.find("debug");
	if(pipe != configMap.end()){
		debug = std::atoi(configMap["debug"].c_str());
		strDebug = configMap["debug"];
	}
	pipe = configMap.find("outputFile");
	if(pipe != configMap.end())
		outputFile = configMap["outputFile"].c_str();
	
	ut = utils(strDebug, outputFile);
	
	pipe = configMap.find("dimensions");
	if(pipe != configMap.end()){
		dim = std::atoi(configMap["dimensions"].c_str());
	}
        configured = true;	
	ut.writeDebug("enumerateUVSCPipe","Configured with parameters { dim: " + std::to_string(dim) + " , debug: " + strDebug + ", outputFile: " + outputFile);
	
	return true;
}

void enumerateUVSCPipe::outputData(pipePacket &inData){
	return;		
}
