#include <string>
#include <vector>
#include <set>
#include <math.h>
#include <algorithm>
#include "alphaComplex.hpp"
#include <fstream>
template<typename nodeType>
alphaComplex<nodeType>::alphaComplex(double maxE, double maxD) : simplexArrayList<nodeType>::simplexArrayList(0,0) {
	this->simplexType = "alphaComplex";
	this->maxEpsilon = maxE;
	this->maxDimension = maxD;
}

template<typename nodeType>
alphaComplex<nodeType>::~alphaComplex(){
	this->simplexList.clear();
}


template<typename nodeType>
std::vector<std::shared_ptr<nodeType>> alphaComplex<nodeType>::getAllDelaunayCofacets(std::shared_ptr<nodeType> simp){
	
	std::vector<std::shared_ptr<nodeType>> ret;
	unsigned dimension  = simp->simplex.size();
	for(auto iter = this->simplexList[dimension].rbegin();iter!=this->simplexList[dimension].rend();++iter){

/*	
	for(auto simplex : this->simplexList[dimension]){
*/		auto simplex = *iter;
		std::vector<unsigned> :: iterator it;
		std::vector<unsigned> v(simplex->simplex.size());
		
		it = std::set_intersection(simp->simplex.begin(),simp->simplex.end(),simplex->simplex.begin(),simplex->simplex.end(),v.begin());
		v.resize(it-v.begin());
		
		if(v.size() == simp->simplex.size()){
			ret.push_back(simplex);
			if(simplex->weight==simp->weight)
				return ret;
		}
	}

	return ret;
}

template<>
void alphaComplex<alphaNode>::buildBetaComplex(std::vector<std::vector<unsigned>> dsimplexmesh, int npts, std::vector<std::vector<double>> inputData,double beta,std::string betaMode){
	unsigned maxDimension = (*dsimplexmesh.begin()).size();
	this->bin = binomialTable(npts, this->maxDimension+1);
	for(int i=0; i <= this->maxDimension; i++)
		this->simplexList.push_back({});
   
	for(auto simplex : dsimplexmesh){
		unsigned int pow_set_size = pow(2, simplex.size());
		for(int counter =1;counter<pow_set_size;counter++){
			double weight =0;
			std::set<unsigned> gensimp;
			for(int j=0;j<simplex.size();j++){
				if(counter & (1<<j)){
					unsigned indnew;
					indnew = *(std::next(simplex.begin(),j));
					for(auto x:gensimp){
						if(x<indnew){
							if(weight<(*(this->distMatrix))[x][indnew])
								weight = (*(this->distMatrix))[x][indnew];
						}
						else if(weight<(*(this->distMatrix))[indnew][x])
							weight = (*(this->distMatrix))[indnew][x];
					}
					gensimp.insert(indnew);
				}
			}

			std::shared_ptr<alphaNode> tot = std::make_shared<alphaNode>(alphaNode(gensimp,weight));
			if(gensimp.size()==1)
				tot->hash = *(gensimp.begin());
			        
			else
				tot->hash = this->simplexHash(gensimp);

			this->simplexList[gensimp.size()-1].insert(tot);
			gensimp.clear();
		}
	}
	
/*
std::vector<std::vector<unsigned>> edges(inputData.size(),std::vector<unsigned>(inputData.size(),0));
for(auto x : this->simplexList[1]){
	std::vector<unsigned> edge;
	for(auto y : x->simplex)
		edge.push_back(y);
    edges[edge[0]][edge[1]] = 1;
	}
this->simplexList.push_back({});

for(auto it = simplexList[maxDimension].begin(); it != simplexList[maxDimension].end(); it++){
		std::set<unsigned> vertices;

		vertices = (*it)->simplex;

		//Iterate over points to possibly add to the simplex
		//Use points larger than the maximal vertex in the simplex to prevent double counting
		unsigned minPt = *vertices.rbegin() + 1;

		for(unsigned pt = minPt; pt < this->simplexList[0].size(); pt++){
			//Compute the weight using all edges
			double maxWeight = (*it)->weight;
			for(auto i : vertices) maxWeight = std::max(maxWeight, (*this->distMatrix)[i][pt]);
			
//i***************************For beta complex valid simplex Condition ****************************			
            bool valid = true;
            if(this->complexType == "alphaComplex"){
                 for(auto i : vertices) 
						if(edges[i][pt]==0){
						   valid = false;
						   break;
					   }
					  
            if(valid){
				std::shared_ptr<alphaNode> tot = std::make_shared<alphaNode>(alphaNode());
				tot->simplex = vertices;
				tot->simplex.insert(pt);
				tot->weight = maxWeight;
				tot->hash = this->simplexHash(tot->simplex);
				simplexList[tot->simplex.size()-1].insert(tot);
				}
//i************************************************************************************************
			}
		}
	}

int di=0;
std::ofstream out("incedenceMatrixGeneraldsphere2"+std::to_string(inputData.size())+betaMode+"beta"+std::to_string(beta)+".csv");

for( auto x : this->simplexList)
	std::cout<<"Count of "<<di++<<"-simplex ::"<<x.size()<<"\n";


for (auto row : edges) {
  for (auto col : row)
    out << col <<' ';
  out << '\n';
}

out.close();

 
	for(auto simplex1 :simplexList[maxDimension+1] ){
	 	auto simplex=simplex1->simplex;
		unsigned int pow_set_size = pow(2, simplex.size());
		for(int counter =1;counter<pow_set_size;counter++){
			double weight =0;
			std::set<unsigned> gensimp;
			for(int j=0;j<simplex.size();j++){
				if(counter & (1<<j)){
					unsigned indnew;
					indnew = *(std::next(simplex.begin(),j));
					for(auto x:gensimp){
						if(x<indnew){
							if(weight<(*(this->distMatrix))[x][indnew])
								weight = (*(this->distMatrix))[x][indnew];
						}
						else if(weight<(*(this->distMatrix))[indnew][x])
							weight = (*(this->distMatrix))[indnew][x];
					}
					gensimp.insert(indnew);
				}
			}

			std::shared_ptr<alphaNode> tot = std::make_shared<alphaNode>(alphaNode(gensimp,weight));
			if(gensimp.size()==1)
				tot->hash = *(gensimp.begin());
			else
				tot->hash = this->simplexHash(gensimp);
			this->simplexList[gensimp.size()-1].insert(tot);
			gensimp.clear();
		}
	}
*/
int di=0;
for( auto x : this->simplexList)
	std::cout<<"Count of "<<di++<<"-simplex ::"<<x.size()<<"\n";
return;
}



template<>
void alphaComplex<alphaNode>::buildBetaComplexFilteration(std::vector<std::vector<unsigned>> dsimplexmesh, int npts, std::vector<std::vector<double>> inputData,kdTree tree){
	unsigned maxDimension = dsimplexmesh[0].size()-1;
	this->bin = binomialTable(npts, this->maxDimension+1);
	for(int i=0; i <= this->maxDimension; i++)
		this->simplexList.push_back({});

     std::cout<<dsimplexmesh.size();
	for(auto simplex : dsimplexmesh){
		unsigned int pow_set_size = pow(2, simplex.size());
		for(int counter =1;counter<pow_set_size;counter++){
			double weight =0;
			std::set<unsigned> gensimp;
			for(int j=0;j<simplex.size();j++){
				if(counter & (1<<j)){
					unsigned indnew;
					indnew = *(std::next(simplex.begin(),j));
					for(auto x:gensimp){
						if(x<indnew){
							if(weight<(*(this->distMatrix))[x][indnew])
								weight = (*(this->distMatrix))[x][indnew];
						}
						else if(weight<(*(this->distMatrix))[indnew][x])
							weight = (*(this->distMatrix))[indnew][x];
					}
					gensimp.insert(indnew);
				}
			}
			//overwriting for Alpha
			double weight1;
			if(gensimp.size()>1)
	      			weight1 = utils::circumRadius(gensimp,this->distMatrix);
   			else
		        	weight1 = weight/2;

			std::shared_ptr<alphaNode> tot = std::make_shared<alphaNode>(alphaNode(gensimp,weight1));
		
			if(gensimp.size()>1)
				tot->circumRadius = utils::circumRadius(gensimp,this->distMatrix);
			else{
		    		tot->circumRadius = weight/2;
			}
			if(gensimp.size()>2)
				tot->circumCenter = utils::circumCenter(gensimp,inputData);
			else if(gensimp.size()==2){
 				auto first = gensimp.begin();
				std::vector<double> R;
				std::vector<double> A = inputData[*first];
      				std::advance(first, 1);
				std::vector<double> B = inputData[*first];
	   			std::transform(A.begin(), A.end(), B.begin(), std::back_inserter(R),[](double e1,double e2){return ((e1+e2)/2);});
				tot->circumCenter = R;
   			 }
			else
   				tot->circumCenter = inputData[*(gensimp.begin())];
			if(gensimp.size()==1)
				tot->hash = *(gensimp.begin());
			else
				tot->hash = this->simplexHash(gensimp);
			this->simplexList[gensimp.size()-1].insert(tot);
			gensimp.clear();
		}
	}

	//ALPHA COMPLEX FILTERTION BASED on FOLLOWING algorithm
	/*Filtration value computation algorithm
	for i : dimension →0 do
	   for all σ of dimension i
	        if filtration(σ) is NaN then
	            filtration(σ)=α2(σ)
	        end if
        	for all τ face of σ do            // propagate alpha filtration value
        	  if  filtration(τ) is not NaN then
        	       filtration(τ) = min( filtration(τ), filtration(σ) )
        	  else
        	     if τ is not Gabriel for σ then
        	       filtration(τ) = filtration(σ)
        	  end if
       	    end if
  	   end for
 	 end for
	end for
	make_filtration_non_decreasing()
	prune_above_filtration()
        */

	for(int dim = this->simplexList.size()-1; dim >= 0; dim--){
		for(auto simplexiter = this->simplexList[dim].rbegin(); simplexiter != this->simplexList[dim].rend(); simplexiter++){
			std::shared_ptr<alphaNode> simplex = (*simplexiter);
			if(simplex->filterationvalue ==-1)
				simplex->filterationvalue = simplex->circumRadius;
			if(dim>0){
				for(auto face : this->simplexList[dim-1]){
					bool gabriel = true;
   			                std::vector<int> neg1;
					std::vector<unsigned> points_check(simplex->simplex.size());
					std::vector<unsigned> guilty_points_check;
					std::vector<unsigned> :: iterator it;
					std::vector<unsigned> v(face->simplex.size());
					it = std::set_intersection(simplex->simplex.begin(),simplex->simplex.end(),face->simplex.begin(),face->simplex.end(),v.begin());
					v.resize(it-v.begin());
					if(v.size() == face->simplex.size()){
						if(face->filterationvalue !=-1)
							face->filterationvalue = std::min(face->filterationvalue, simplex->filterationvalue);				
						else {

			                                std::vector<size_t> neighborsface = tree.neighborhoodIndices(face->circumCenter, face->circumRadius); //All neighbors in epsilon-ball
 					   	        std::vector<unsigned> toremove(face->simplex.begin(),face->simplex.end()); 
					   	        std::sort(toremove.begin(),toremove.end());
    						        std::sort(neighborsface.begin(), neighborsface.end());
   						        std::set_difference(neighborsface.begin(),neighborsface.end(),toremove.begin(),toremove.end(), std::back_inserter(neg1));
							/*
							std::vector<unsigned>::iterator it;											
							it=std::set_difference (simplex->simplex.begin(), simplex->simplex.end(), face->simplex.begin(), face->simplex.end(), points_check.begin());
							points_check.resize(it-points_check.begin());
							for(it = points_check.begin(); it !=points_check.end();++it){
								std::vector<double> coordinates;
								for(int i =0;i<face->circumCenter.size();i++)
									coordinates.push_back(inputData[*it][i]);						
								double distance = utils::vectors_distance(coordinates,face->circumCenter);
								if(pow(distance,2)<face->circumRadius)	
									gabriel = false;
								guilty_points_check.push_back((*it));
							}
							*/
							if(neg1.size()!=0)
								gabriel = false;
						}   
					}
					if(!gabriel){
						std::vector<unsigned> v(neg1.size()+face->simplex.size());
						std::vector<unsigned>::iterator it;
						it=std::set_union (face->simplex.begin(), face->simplex.end(), neg1.begin(), neg1.end(), v.begin());
						v.resize(it-v.begin());
						std::sort(v.begin(),v.end());
					/*	for(auto face1 : this->simplexList[v.size()-1]){
							std::vector<unsigned> v1(simplex->simplex.size());
							it = std::set_intersection(v.begin(), v.end(), face1->simplex.begin(), face1->simplex.end(), v1.begin());
							v1.resize(it-v1.begin());
							if(v1.size() == face1->simplex.size())
							*/
                                                std::set<unsigned> simpnode(v.begin(),v.end());
						if(simpnode.size() <= this->simplexList.size()){
						if(simpnode.size()>1)
							face->filterationvalue = utils::circumRadius(simpnode,this->distMatrix);
						else
					        	face->filterationvalue = 0;
						}
						else
						{
								
		        int n = simpnode.size();
			int dim = this->simplexList.size();
			std::vector<unsigned> dsimplex(dim);
			std::set<unsigned> dsimplexIndexed;
			double radius;
                        double maxradius = 0;
			std::vector<unsigned>::iterator first = dsimplex.begin(), last = dsimplex.end();

			std::generate(first, last, UniqueNumber);
			
			
			for(int i=0;i<dim;i++)
				dsimplexIndexed.insert(v[dsimplex[i]-1]);
                       
			radius = utils::circumRadius(dsimplexIndexed,this->distMatrix);
			if(radius>maxradius)
			       maxradius=radius;	
            
	
	while((*first) != n-dim+1){  	

  	    std::vector<unsigned>::iterator mt = last;
        while (*(--mt) == n-(last-mt)+1);
        (*mt)++;
        while (++mt != last) *mt = *(mt-1)+1;

        std::set<unsigned> dsimplexIndexed1;
       for(int i=0;i<dim;i++)
			dsimplexIndexed1.insert(v[dsimplex[i]-1]);
            
			radius = utils::circumRadius(dsimplexIndexed1,this->distMatrix);
			if(radius>maxradius)
			       maxradius=radius;	
    	 	}      
	       face->filterationvalue = maxradius;
	  
						}
					}
				}
			}
		}
	}
	
	//Reinserting to sort by filterationvalue and remove simplexes with weight greater than alphafilteration value (maxEpsilon)
	std::vector<std::set<std::shared_ptr<alphaNode>, cmpByWeight<std::shared_ptr<alphaNode>>>> simplexList1;		//Holds ordered list of simplices in each dimension
	for(int dim=0;dim < this->simplexList.size();dim++){
		simplexList1.push_back({});
	   	for(auto simplex : this->simplexList[dim]){
			 if(simplex->filterationvalue <= this->alphaFilterationValue){ //Valid Simplex after filteration
				 simplex->weight = simplex->filterationvalue;
				 simplexList1[dim].insert(simplex);
			 }
		}
	}
	this->simplexList = simplexList1;
	return;
}



template<>
void alphaComplex<alphaNode>::buildAlphaComplex(std::vector<std::vector<unsigned>> dsimplexmesh, int npts, std::vector<std::vector<double>> inputData){
	unsigned maxDimension = dsimplexmesh[0].size()-1;
	this->bin = binomialTable(npts, this->maxDimension+1);
	for(int i=0; i <= this->maxDimension; i++)
		this->simplexList.push_back({});

     std::cout<<dsimplexmesh.size();
	for(auto simplex : dsimplexmesh){
		unsigned int pow_set_size = pow(2, simplex.size());
		for(int counter =1;counter<pow_set_size;counter++){
			double weight =0;
			std::set<unsigned> gensimp;
			for(int j=0;j<simplex.size();j++){
				if(counter & (1<<j)){
					unsigned indnew;
					indnew = *(std::next(simplex.begin(),j));
					for(auto x:gensimp){
						if(x<indnew){
							if(weight<(*(this->distMatrix))[x][indnew])
								weight = (*(this->distMatrix))[x][indnew];
						}
						else if(weight<(*(this->distMatrix))[indnew][x])
							weight = (*(this->distMatrix))[indnew][x];
					}
					gensimp.insert(indnew);
				}
			}
			//overwriting for Alpha
			double weight1;
			if(gensimp.size()>1)
	      			weight1 = utils::circumRadius(gensimp,this->distMatrix);
   			else
		        	weight1 = weight/2;

			std::shared_ptr<alphaNode> tot = std::make_shared<alphaNode>(alphaNode(gensimp,weight1));
		
			if(gensimp.size()>1)
				tot->circumRadius = utils::circumRadius(gensimp,this->distMatrix);
			else{
		    		tot->circumRadius = weight/2;
			}
			if(gensimp.size()>2)
				tot->circumCenter = utils::circumCenter(gensimp,inputData);
			else if(gensimp.size()==2){
 				auto first = gensimp.begin();
				std::vector<double> R;
				std::vector<double> A = inputData[*first];
      				std::advance(first, 1);
				std::vector<double> B = inputData[*first];
	   			std::transform(A.begin(), A.end(), B.begin(), std::back_inserter(R),[](double e1,double e2){return ((e1+e2)/2);});
				tot->circumCenter = R;
   			 }
			else
   				tot->circumCenter = inputData[*(gensimp.begin())];
			if(gensimp.size()==1)
				tot->hash = *(gensimp.begin());
			else
				tot->hash = this->simplexHash(gensimp);
			this->simplexList[gensimp.size()-1].insert(tot);
			gensimp.clear();
		}
	}

	//ALPHA COMPLEX FILTERTION BASED on FOLLOWING algorithm
	/*Filtration value computation algorithm
	for i : dimension →0 do
	   for all σ of dimension i
	        if filtration(σ) is NaN then
	            filtration(σ)=α2(σ)
	        end if
        	for all τ face of σ do            // propagate alpha filtration value
        	  if  filtration(τ) is not NaN then
        	       filtration(τ) = min( filtration(τ), filtration(σ) )
        	  else
        	     if τ is not Gabriel for σ then
        	       filtration(τ) = filtration(σ)
        	  end if
       	    end if
  	   end for
 	 end for
	end for
	make_filtration_non_decreasing()
	prune_above_filtration()
        */

	for(int dim = this->simplexList.size()-1; dim >= 0; dim--){
		for(auto simplexiter = this->simplexList[dim].rbegin(); simplexiter != this->simplexList[dim].rend(); simplexiter++){
			std::shared_ptr<alphaNode> simplex = (*simplexiter);
			if(simplex->filterationvalue ==-1)
				simplex->filterationvalue = simplex->circumRadius;
			if(dim>0){
				for(auto face : this->simplexList[dim-1]){
					bool gabriel = true;
					std::vector<unsigned> points_check(simplex->simplex.size());
					std::vector<unsigned> guilty_points_check;
					std::vector<unsigned> :: iterator it;
					std::vector<unsigned> v(face->simplex.size());
					it = std::set_intersection(simplex->simplex.begin(),simplex->simplex.end(),face->simplex.begin(),face->simplex.end(),v.begin());
					v.resize(it-v.begin());
					if(v.size() == face->simplex.size()){
						if(face->filterationvalue !=-1)
							face->filterationvalue = std::min(face->filterationvalue, simplex->filterationvalue);				
						else {
							std::vector<unsigned>::iterator it;											
							it=std::set_difference (simplex->simplex.begin(), simplex->simplex.end(), face->simplex.begin(), face->simplex.end(), points_check.begin());
							points_check.resize(it-points_check.begin());
							for(it = points_check.begin(); it !=points_check.end();++it){
								std::vector<double> coordinates;
								for(int i =0;i<face->circumCenter.size();i++)
									coordinates.push_back(inputData[*it][i]);						
								double distance = utils::vectors_distance(coordinates,face->circumCenter);
								if(pow(distance,2)<face->circumRadius)	
									gabriel = false;
								guilty_points_check.push_back((*it));
							}
						}
					}
					if(!gabriel){
						std::vector<unsigned> v(simplex->simplex.size());
						std::vector<unsigned>::iterator it;
						it=std::set_union (face->simplex.begin(), face->simplex.end(), guilty_points_check.begin(), guilty_points_check.end(), v.begin());
						v.resize(it-v.begin());
						for(auto face1 : this->simplexList[v.size()-1]){
							std::vector<unsigned> v1(simplex->simplex.size());
							it = std::set_intersection(v.begin(), v.end(), face1->simplex.begin(), face1->simplex.end(), v1.begin());
							v1.resize(it-v1.begin());
							if(v1.size() == face1->simplex.size())
								face->filterationvalue = face1->filterationvalue;
						}
					}
				}
			}
		}
	}
	
	//Reinserting to sort by filterationvalue and remove simplexes with weight greater than alphafilteration value (maxEpsilon)
	std::vector<std::set<std::shared_ptr<alphaNode>, cmpByWeight<std::shared_ptr<alphaNode>>>> simplexList1;		//Holds ordered list of simplices in each dimension
	for(int dim=0;dim < this->simplexList.size();dim++){
		simplexList1.push_back({});
	   	for(auto simplex : this->simplexList[dim]){
			 if(simplex->filterationvalue <= this->alphaFilterationValue){ //Valid Simplex after filteration
				 simplex->weight = simplex->filterationvalue;
				 simplexList1[dim].insert(simplex);
			 }
		}
	}
	this->simplexList = simplexList1;
	return;
}




template<typename nodeType>
void alphaComplex<nodeType>::buildBetaComplexFilteration(std::vector<std::vector<unsigned>> dsimplexmesh, int npts, std::vector<std::vector<double>> inputData,kdTree tree){	
	this->ut.writeLog(this->simplexType,"No build alpha complex function defined");
	return;
}

template<typename nodeType>
void alphaComplex<nodeType>::buildAlphaComplex(std::vector<std::vector<unsigned>> dsimplexmesh, int npts, std::vector<std::vector<double>> inputData){	
	this->ut.writeLog(this->simplexType,"No build alpha complex function defined");
	return;
}



template<typename nodeType>
void alphaComplex<nodeType>::buildBetaComplex(std::vector<std::vector<unsigned>> dsimplexmesh, int npts, std::vector<std::vector<double>> inputData,double beta,std::string betaMode){
	this->ut.writeLog(this->simplexType,"No build beta complex function defined");
	return;
}


//Explicit Template Class Instantiation
template class alphaComplex<simplexNode>;
template class alphaComplex<alphaNode>;
template class alphaComplex<witnessNode>;

