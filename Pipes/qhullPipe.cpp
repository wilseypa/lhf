
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
#include "alphaComplex.hpp"
#include "betaComplex.hpp"


template <typename nodeType>
qhullPipe<nodeType>::qhullPipe(){
    /**
	    qhullPipe()
	 
		@brief Class constructor
		@tparam nodeType The data type of the simplex node.
	*/
	this->pipeType = "qhullPipe";
	return;
}


template <typename nodeType>
void qhullPipe<nodeType>::runPipe(pipePacket<nodeType> &inData){
	/**
		runPipe(pipePacket<nodeType> &inData)
		
		@brief Run qhull pipe to generate delaunay-based complexes
		@tparam nodeType The data type of the simplex node.
		@param inData The pipepacket data holding working data for qhull.
	*/
	
    Qhull qh;
    
    //TODO: Is this kdtree below required for all 3 of these complexes? I assume not...
    //kdTree tree(inData.inputData, inData.inputData.size()); //KDTree for efficient nearest neighbor search
   
    std::vector<double> sdata;
    
    //serializing all the data
    for(auto a : inData.inputData)
 			for(auto b : a)
 				sdata.push_back(b);
 	
 	//PointCoordinates data type defined in ...
    PointCoordinates *pts = new PointCoordinates(qh,inData.inputData[0].size(),"UCI Data Sets");
    pts->append(sdata);
    
    
    if(this->mode == "alpha" || this->mode == "beta"){
    
        qh.runQhull(pts->comment().c_str(),pts->dimension(),pts->count(),&*pts->coordinates(),"d o");
        qdelaunay_o(qh, inData.complex->dsimplexmesh);    
            
        if(this->mode == "beta"){
            //Beta complex uses sparsification filter for B < 1 and B > 1
            //((betaComplex<nodeType>*)inData.complex)->buildBetaComplexFilteration(dsimplexes, inData.inputData.size(),inData.inputData, tree);
            //((betaComplex<nodeType>*)inData.complex)->buildBetaComplex(dsimplexes, inData.inputData.size(),inData.inputData,1,"highDim");
            std::cout << "TODO: Beta-sparsified Complex" << std::endl;
        }else { //alpha
            //Alpha complex uses gabriel filtration (circumradius at most alpha)
            ((alphaComplex<nodeType>*)inData.complex)->buildAlphaComplex(inData.complex->dsimplexmesh,inData.inputData.size(),inData.inputData);
        }
    } else if(this->mode == "weightedAlpha"){
		qh.runQhull(pts->comment().c_str(),pts->dimension(),pts->count(),&*pts->coordinates(),"d Qt");
        qdelaunay_o(qh, inData.complex->dsimplexmesh);    
        
        ((alphaComplex<nodeType>*)inData.complex)->buildWeightedAlphaComplex(inData.complex->dsimplexmesh,inData.inputData.size(),inData.inputData);
	}

    //If we are not in debug mode, clear the simplex mesh here; in debug mode, output and then clear.
    if(this->debug == 0)
        inData.complex->dsimplexmesh = std::vector<std::vector<unsigned>>();

	this->ut.writeDebug("qhullPipe", "\tSuccessfully Executed pipe");
	return;
}

template <typename nodeType>
void qhullPipe<nodeType>::qdelaunay_o(const Qhull &qhull, std::vector<std::vector<unsigned>> &regions){
    /**
		qdelaunay_o(const Qhull &qhull, std::vector<std::vector<unsigned>> &regions)
		
		@brief Wrapper to call qhull and store simplex mesh into the complex
		@tparam nodeType The data type of the simplex node.
		@param qhull The qhull class after executing runQhull
        @param regions Pointer to the complex simplexMesh for output
	*/
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

	QhullFacetListIterator k(facets);
	while(k.hasNext()){
		QhullFacet f = k.next();
		std::vector<unsigned> vertices;
		if(!f.isUpperDelaunay()){
			if(!f.isTopOrient() && f.isSimplicial()){
				QhullVertexSet vs = f.vertices();
				for(int i=0;i<(int)vs.size();++i){
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
    return;
}


// configPipe -> configure the function settings of this pipeline segment
template <typename nodeType>
bool qhullPipe<nodeType>::configPipe(std::map<std::string, std::string> &configMap){
    /**
	    configPipe(std::map<std::string, std::string> &configMap)
	 
		@brief Configures the pipe and sets arguments based on the configMap passed. Called before execution (runPipe). If required values not found or configuration is invalid, returns false. 
		@tparam nodeType The data type of the simplex node.
		@param configMap The configuration map for this pipeline
        @return boolean
	*/
	std::string strDebug;

	auto pipe = configMap.find("debug");
	if(pipe != configMap.end()){
		this->debug = std::atoi(configMap["debug"].c_str());
		strDebug = configMap["debug"];
	}
	pipe = configMap.find("outputFile");
	if(pipe != configMap.end())
		this->outputFile = configMap["outputFile"].c_str();
		
	pipe = configMap.find("mode");
	if(pipe != configMap.end())
		this->mode = configMap["mode"].c_str();

	this->ut = utils(strDebug, this->outputFile);

	this->configured = true;
	this->ut.writeDebug("qhullPipe","Configured with parameters { eps: " + configMap["epsilon"] + " , debug: " + strDebug + ", outputFile: " + this->outputFile + ", Mode: " + this->mode + " }");

	return true;
}

// outputData -> used for tracking each stage of the pipeline's data output without runtime
template <typename nodeType>
void qhullPipe<nodeType>::outputData(pipePacket<nodeType> &inData){
    /**
	    outputData(pipePacket<nodeType> &inData)
	 
		@brief Outputs dsimplexmesh to a file if debug mode is true. 
		@tparam nodeType The data type of the simplex node.
		@param inData The pipePacket data being used in the pipeline.
	*/
	std::ofstream file;
	file.open("output/" + this->pipeType + "_output.csv");
        
    for(auto a : inData.complex->dsimplexmesh){
		for(auto d : a){
			file << d << ",";
		}
		file << "\n";
	}

	file.close();
    
    //Clear the dsimplexmesh
    inData.complex->dsimplexmesh = std::vector<std::vector<unsigned>>();
    
	file.open("output/" + inData.complex->complexType + "_output.csv");
        
    for(auto a : inData.complex->simplexList){
		for(auto d : a){
			file << d->weight << ",[ ";
            for (auto ind : d->simplex){
                file << ind << " ";
            }
            file << "]\n";
		}
	}

	file.close();
    
	return;
}

template class qhullPipe<simplexNode>;
template class qhullPipe<alphaNode>;
template class qhullPipe<witnessNode>;
