
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
#include "qhullPipe.hpp"
#include "utils.hpp"

// basePipe constructor
qhullPipe::qhullPipe(){
	pipeType = "qhullPipe";
	return;
}

// runPipe -> Run the configured functions of this pipeline segment
void qhullPipe::runPipe(pipePacket &inData){
    Qhull qh;
    std::vector<double> sdata;
    //serializing all the data
    for(auto a : inData.inputData)
 			for(auto b : a)
 				sdata.push_back(b);
    PointCoordinates *pts = new PointCoordinates(qh,inData.inputData[0].size(),"UCI Data Sets");
    pts->append(sdata);
    qh.runQhull(pts->comment().c_str(),pts->dimension(),pts->count(),&*pts->coordinates(),"d o");
    auto dsimplexes = qdelaunay_o(qh);
    inData.complex->buildAlphaComplex(dsimplexes,inData.inputData.size(),inData.inputData);

	ut.writeDebug("qhullPipe", "\tSuccessfully Executed pipe");
	return;
}

std::vector<std::vector<int>>  qhullPipe::qdelaunay_o(const Qhull &qhull){
	int hullDimension = qhull.hullDimension();
        std::vector<std::vector<double> > inputSites;
	QhullPoints points = qhull.points();

	QhullPointsIterator j(points);
	while(j.hasNext()){
		QhullPoint point = j.next();
		inputSites.push_back(point.toStdVector());
	}
	QhullFacetList facets = qhull.facetList();
	int numFacets = facets.count();
	size_t numRidges = numFacets*hullDimension/2;

	std::vector<std::vector<int>> regions;
	QhullFacetListIterator k(facets);
	while(k.hasNext()){
		QhullFacet f = k.next();
		std::vector<int> vertices;
		if(!f.isUpperDelaunay()){
			if(!f.isTopOrient() && f.isSimplicial()){
				QhullVertexSet vs = f.vertices();
				vertices.push_back(vs[1].point().id());
				vertices.push_back(vs[0].point().id());
				for(int i=2;i<(int)vs.size();++i){
					vertices.push_back(vs[i].point().id());
				}
			}
			else{
				QhullVertexSetIterator i(f.vertices());
				while(i.hasNext()){
					QhullVertex vertex = i.next();
					QhullPoint p = vertex.point();
					vertices.push_back(p.id());
				}
			}
			regions.push_back(vertices);
		}
	}
      return regions;
}


// configPipe -> configure the function settings of this pipeline segment
bool qhullPipe::configPipe(std::map<std::string, std::string> &configMap){
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

	configured = true;
	ut.writeDebug("qhullPipe","Configured with parameters { eps: " + configMap["epsilon"] + " , debug: " + strDebug + ", outputFile: " + outputFile + " }");

	return true;
}

// outputData -> used for tracking each stage of the pipeline's data output without runtime
void qhullPipe::outputData(pipePacket &inData){
	std::ofstream file;
	file.open("output/" + pipeType + "_output.csv");



	file.close();
	return;
}
