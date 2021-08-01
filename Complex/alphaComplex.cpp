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
std::vector<std::shared_ptr<nodeType>> alphaComplex<nodeType>::getAllDelaunayCofacets(std::shared_ptr<nodeType> simp, std::unordered_map<std::shared_ptr<nodeType>,std::shared_ptr<nodeType>> pivotPairs,bool emergent){	//For each pivot, which column has that pivot
	
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
			if(emergent&&simplex->weight==simp->weight){
				if(pivotPairs.find(simplex) == pivotPairs.end()) return ret;
					emergent = false;
			}
		}
	}

	return ret;
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
		
		if(v.size() == simp->simplex.size())
			ret.push_back(simplex);
	}

	return ret;
}
template<typename nodeType>

bool  alphaComplex<nodeType>::checkGabriel(std::vector<double> point, std::vector<unsigned> dsimplex ,std::vector<std::vector<double>> &inputData, double beta)
{
	bool intersectionCircle= false;
	bool intersectionLune = false;
	if(beta <0)
		exit(0);
	else if(beta==0)
		return false;
	else if(beta <1)
		intersectionCircle = true;
	else
		intersectionLune = true;

	auto betacentersandradii = utils::calculateBetaCentersandRadius(dsimplex ,inputData,this->distMatrix,beta);
	int i=0;
    	for(auto bc :betacentersandradii.first){
			double distance = utils::vectors_distance(point,bc);
			if(intersectionCircle&&distance > betacentersandradii.second[i])
				return false;
			if(!intersectionCircle&&distance < betacentersandradii.second[i])
				return true;
       i++;
		}
	if(intersectionCircle)
		return true;
	else
		return false;
}
template<>
void alphaComplex<alphaNode>::buildFilteration(std::vector<std::vector<unsigned>> dsimplexmesh, int npts, std::vector<std::vector<double>> inputData,double beta,kdTree tree){
	unsigned maxDimension = dsimplexmesh[0].size()-1;
	this->bin = binomialTable(npts, this->maxDimension+1);
	for(int i=0; i <= maxDimension; i++)
		this->simplexList.push_back({});

int pk;
	for(auto simplex : dsimplexmesh){
		unsigned int pow_set_size = pow(2, simplex.size());
		for(int counter =1;counter<pow_set_size;counter++){
			double weight =0;
			std::set<unsigned> gensimp;
			for(int j=0;j<simplex.size();j++){
				if(counter & (1<<j)){
					unsigned indnew;
					indnew = *(std::next(simplex.begin(),j));
					gensimp.insert(indnew);
				}
			}
			std::shared_ptr<alphaNode> tot = std::make_shared<alphaNode>(alphaNode(gensimp,0));
			if(gensimp.size()==1)
				tot->hash = *(gensimp.begin());
			        
			else
				tot->hash = this->simplexHash(gensimp);

			if(gensimp.size()>1)
				tot->circumRadius = utils::circumRadius(gensimp,this->distMatrix);
			else{
		    		tot->circumRadius = pow(weight/2,2);
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

			this->simplexList[gensimp.size()-1].insert(tot);
			gensimp.clear();
		}
	}

	for(auto dim = maxDimension;dim >0;dim--)
		for(auto simplex : this->simplexList[dim]){
			std::vector<unsigned> vecsimplex(simplex->simplex.begin(),simplex->simplex.end());
			if(simplex->weight == 0)
				simplex->weight = (dim>1? utils::circumRadius(simplex->simplex,this->distMatrix):(dim==1?(pow((*(this->distMatrix))[std::min(vecsimplex[0],vecsimplex[1])][std::max(vecsimplex[0],vecsimplex[1])],2)):0));
                 
                 
         	std::vector<std::set<std::shared_ptr<alphaNode>, cmpByWeight<std::shared_ptr<alphaNode>>>> simplexfaceList;
          	 for(int i=0;i<simplex->simplex.size();i++)
			  simplexfaceList.push_back({}); 
		unsigned int pow_set_size = pow(2, simplex->simplex.size());
		for(int counter =1;counter<pow_set_size;counter++){
			double weight =0;
			std::set<unsigned> gensimp;
			for(int j=0;j<simplex->simplex.size();j++){
				if(counter & (1<<j)){
					unsigned indnew;
					indnew = *(std::next(vecsimplex.begin(),j));
					gensimp.insert(indnew);
				}
			}
				for(auto x:this->simplexList[gensimp.size()-1])
					if(x->simplex==gensimp){
						simplexfaceList[gensimp.size()-1].insert(x);
						break;
		       			}
		
		}
                for(int facedim = simplex->simplex.size()-2;facedim>0;facedim--){
			for(auto face : simplexfaceList[facedim])
				if(face->weight!=0){
					for(auto coface:simplexList[facedim+1]){
						std::vector<unsigned> A(coface->simplex.begin(),coface->simplex.end());
						std::vector<unsigned> B(face->simplex.begin(),face->simplex.end());
						if(utils::isSubset(A,B)){
						   face->weight = std::min(coface->weight,face->weight);
						}
					}
				}
				else
				{
					for(auto coface:simplexList[facedim+1]){
						std::vector<unsigned> A(coface->simplex.begin(),coface->simplex.end());
						std::vector<unsigned> B(face->simplex.begin(),face->simplex.end());
						if(utils::isSubset(A,B)){
							std::vector<unsigned> points_check(A.size());	
							std::vector<unsigned>::iterator it;
							it=std::set_difference (A.begin(),A.end(), B.begin(), B.end(), points_check.begin());
							points_check.resize(it-points_check.begin());
                                                        std::vector<double> coordinates;
							for(it = points_check.begin(); it !=points_check.end();++it){
								std::vector<double> coordinates;
								for(int i =0;i<face->circumCenter.size();i++)
									coordinates.push_back(inputData[*it][i]);					                                         
								double distance = utils::vectors_distance(coordinates,face->circumCenter);
						//		if(distance < sqrt(face->circumRadius))
								if(checkGabriel(coordinates,B,inputData,beta))
									if(face->weight ==0)
										face->weight = coface->weight;
									else
										face->weight = std::min(face->weight,coface->weight);
							}
						}
					}
				}
			}
                  for(auto x : simplexfaceList)
			  for(auto simp : x){
				  for(auto relList : this->simplexList[simp->simplex.size()-1])
					  if(simp->simplex == relList->simplex)
						  if(relList ->weight > simp->weight || relList->weight==0)
							  if(simp->simplex.size()>1)
							 	 relList->weight = simp->weight;
			  }
				  


		}

//Simplex List after gabrel weight assignment; now make it non-decreasing

	for(auto x:this->simplexList)
	if(!((*x.begin())->simplex.size()>=this->simplexList.size()))
		for(auto y : x)
		{
			for(auto simplex: simplexList[y->simplex.size()]){	
				std::vector<unsigned> A(simplex->simplex.begin(),simplex->simplex.end());
				std::vector<unsigned> B(y->simplex.begin(),y->simplex.end());
				if(utils::isSubset(A,B)){
					if(simplex->weight < y->weight)
						simplex->weight = y->weight;
				}
			}
		}

	std::vector<std::set<std::shared_ptr<alphaNode>, cmpByWeight<std::shared_ptr<alphaNode>>>> simplexList1;		//Holds ordered list of simplices in each dimension
	for(int dim=0;dim < this->simplexList.size();dim++){
		simplexList1.push_back({});
	   	for(auto simplex : this->simplexList[dim]){
			 if(simplex->weight <= this->alphaFilterationValue){ 
				 simplexList1[dim].insert(simplex);
			 }
		}

	}

	this->simplexList = simplexList1;
int di=0;
for( auto x : this->simplexList)
	std::cout<<"Count of "<<di++<<"-simplex ::"<<x.size()<<"\n";
return;
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
                         
		
			if(gensimp.size()>1)
	      			weight = utils::circumRadius(gensimp,this->distMatrix);
   			else
		        	weight = weight/2;
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
				tot->circumRadius = sqrt(utils::circumRadius(gensimp,this->distMatrix));
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
			if(gensimp.size() == dsimplexmesh[0].size()){
				tot->filterationvalue = tot->circumRadius;
				std::cout<<" "<<tot->circumRadius<<"\n";
			}
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
int cnt=0;
	for(int dim = this->simplexList.size()-1; dim >= 0; dim--){
		for(auto simplexiter = this->simplexList[dim].rbegin(); simplexiter != this->simplexList[dim].rend(); ++simplexiter){
			std::shared_ptr<alphaNode> simplex = (*simplexiter);
			if(simplex->filterationvalue==0)
				simplex->filterationvalue = simplex->circumRadius;
			if(dim>0){for(int d = dim-1;d>=0;d--)
				for(auto face : this->simplexList[d]){
					bool gabriel = true;
					std::vector<unsigned> A(simplex->simplex.begin(),simplex->simplex.end());
					std::vector<unsigned> B(face->simplex.begin(),face->simplex.end());
			                std::vector<size_t> neighborsface;
					std::sort(A.begin(),A.end());
					std::sort(B.begin(),B.end());
                                        if(utils::isSubset(A,B)){
						if(face->filterationvalue !=0)
							face->filterationvalue = std::min(face->filterationvalue, simplex->filterationvalue);				
						else {
							neighborsface.clear();
			                                neighborsface = tree.neighborhoodIndices(face->circumCenter, face->circumRadius); //All neighbors in epsilon-ball
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
							std::vector<size_t> neg1(neighborsface.begin(),neighborsface.end());
							for( auto x : face->simplex)
						         	neg1.erase(std::remove(neg1.begin(), neg1.end(), x),neg1.end());
							neighborsface = neg1;
							if(neg1.size()>0)
								gabriel = false;
						}   
					}	
					if(!gabriel){
						std::vector<double> weighttoselect;
						double weightassigned;
						if(face->filterationvalue !=0)
							if(face->filterationvalue>face->circumRadius)
								weightassigned = face->filterationvalue;
							else
								weightassigned = simplex->filterationvalue;
						else
							weightassigned = simplex->filterationvalue;
						unsigned chosen_coface=-1 ;
						for(unsigned i :neighborsface){
							std::vector<unsigned> k;
							k.push_back(i);
					
						std::vector<unsigned> face1(face->simplex.begin(),face->simplex.end());
							if(!utils::isSubset(face1,k)){
						face1.push_back(i);
							for(auto iter =this->simplexList[face1.size()-1].rbegin();iter != this->simplexList[face1.size()-1].rend();++iter){
                                        	       	        auto cursimplex = (*iter);
								std::vector<unsigned> simplexvertex(cursimplex->simplex.begin(),cursimplex->simplex.end());
								std::sort(simplexvertex.begin(),simplexvertex.end());
								std::sort(face1.begin(),face1.end());
								if(simplexvertex== face1)
								{
									weighttoselect.push_back(cursimplex->filterationvalue);
									 if(weightassigned > cursimplex->filterationvalue){
										 weightassigned = cursimplex->filterationvalue;
										 chosen_coface = i;
										 cursimplex->weight = -1;
									 }
								 }
							}
							}
						}
					/*	for(unsigned i :neighborsface){
						std::vector<unsigned> face1(face->simplex.begin(),face->simplex.end());
						face1.push_back(i);
						std::sort(face1.begin(),face1.end());
							for(auto iter =this->simplexList[face1.size()-1].rbegin();iter != this->simplexList[face1.size()-1].rend();++iter){
                                        	       	        auto cursimplex = (*iter);
								std::vector<unsigned> simplexvertex(cursimplex->simplex.begin(),cursimplex->simplex.end());
								std::sort(simplexvertex.begin(),simplexvertex.end());
								std::sort(face1.begin(),face1.end());
								if(simplexvertex== face1&&i!=chosen_coface){
									int p=0;
									std::cout<<"i=="<<i;
									cursimplex->filterationvalue = -10;
								}
							}
						}
                                       */
						face->filterationvalue = weightassigned;
					/*	auto minmax = std::minmax_element(weighttoselect.begin(),weighttoselect.end()); 
						if(*minmax.first > face->circumRadius)
							face->filterationvalue  =*minmax.first;
						else
							face->filterationvalue =  face->circumRadius;
							*/
					//	face->filterationvalue = minmax.first;
						std::cout<<"\nNOT Gabriel Assigned  weight "<<dim-1<<" "<<face->filterationvalue;
					}
				//	face->filterationvalue=face->circumRadius;
				}
			}
		}
	}
/*
bool done = false;
for(auto x:simplexList){
	for(auto y :x){
		if(y->simplex.size()==simplexList.size()){
			done = true;
			break;
		}
		else
		for(auto simplex : simplexList[y->simplex.size()]){
			std::vector<unsigned> A (simplex->simplex.begin(),simplex->simplex.end());
			std::vector<unsigned> B (y->simplex.begin(),y->simplex.end());
			if(utils::isSubset(A,B)){
                               if(y->filterationvalue==-10) 
			        	simplex->filterationvalue=std::min(simplex->filterationvalue,y->filterationvalue);
			       else if(simplex->weight!=-1)
			        	simplex->filterationvalue=std::max(simplex->filterationvalue,y->filterationvalue);
			       else
				       simplex
		 	}
		}
	}
		if(done)
			break;
	}
	*/
	
	//Reinserting to sort by filterationvalue and remove simplexes with weight greater than alphafilteration value (maxEpsilon)
	std::vector<std::set<std::shared_ptr<alphaNode>, cmpByWeight<std::shared_ptr<alphaNode>>>> simplexList1;		//Holds ordered list of simplices in each dimension
	for(int dim=0;dim < this->simplexList.size();dim++){
		simplexList1.push_back({});
	   	for(auto simplex : this->simplexList[dim]){
			 if(simplex->filterationvalue <= this->alphaFilterationValue){ 
				 simplex->weight = simplex->filterationvalue;
				 simplexList1[dim].insert(simplex);
			 }
		}
	}
	this->simplexList = simplexList1;


int di=0;
for( auto x : this->simplexList)
	std::cout<<"Count of "<<di++<<"-simplex ::"<<x.size()<<"\n";
	return;

}



template<>
void alphaComplex<alphaNode>::buildAlphaComplex(std::vector<std::vector<unsigned>> dsimplexmesh, int npts, std::vector<std::vector<double>> inputData){
	unsigned maxDimension = dsimplexmesh[0].size()-1;
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
			//overwriting for Alpha
			double weight1;
			if(gensimp.size()>1)
	      			weight1 = utils::circumRadius(gensimp,this->distMatrix);
   			else
		        	weight1 = weight/2;

			std::shared_ptr<alphaNode> tot = std::make_shared<alphaNode>(alphaNode(gensimp,weight1));
		
			if(gensimp.size()>1)
				tot->circumRadius = sqrt(utils::circumRadius(gensimp,this->distMatrix));
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
		for(auto simplexiter = this->simplexList[dim].rbegin(); simplexiter != this->simplexList[dim].rend(); ++simplexiter){
			std::shared_ptr<alphaNode> simplex = (*simplexiter);
			if(simplex->filterationvalue ==0)
				simplex->filterationvalue = simplex->circumRadius;
			if(dim>0){for(int d = dim-1;d>=dim-1;d--)
				for(auto face : this->simplexList[d]){
	
					bool gabriel = true;
					std::vector<unsigned> points_check(simplex->simplex.size());
					std::vector<unsigned> guilty_points_check;
					std::vector<unsigned> A(simplex->simplex.begin(),simplex->simplex.end());
					std::vector<unsigned> B(face->simplex.begin(),face->simplex.end());
					std::sort(A.begin(),A.end());
					std::sort(B.begin(),B.end());
					if(utils::isSubset(A,B)){
						if(face->filterationvalue !=0)
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
								if(distance<face->circumRadius){
									gabriel = false;
									guilty_points_check.push_back((*it));
								}
							}
						}
					}
					if(!gabriel){
						std::vector<unsigned> v(simplex->simplex.size());
						std::vector<unsigned>::iterator it;
						
						it=std::set_union (face->simplex.begin(), face->simplex.end(), guilty_points_check.begin(), guilty_points_check.end(), v.begin());
						v.resize(it-v.begin());
						for(auto iter =this->simplexList[v.size()-1].rbegin();iter != this->simplexList[v.size()-1].rend();++iter){
							auto face1 = (*iter);
							std::vector<unsigned> v1(simplex->simplex.size());
							it = std::set_intersection(v.begin(), v.end(), face1->simplex.begin(), face1->simplex.end(), v1.begin());
							v1.resize(it-v1.begin());
							if(v1.size() == face1->simplex.size()){
								face->filterationvalue = face1->filterationvalue;
							}
						}
					}
				
				}
			}
		}
	}



	for(auto x:this->simplexList)
	if(!((*x.begin())->simplex.size()>=this->simplexList.size()))
		for(auto y : x)
		{
			for(auto simplex: simplexList[y->simplex.size()]){	
				std::vector<unsigned> A(simplex->simplex.begin(),simplex->simplex.end());
				std::vector<unsigned> B(y->simplex.begin(),y->simplex.end());
				if(utils::isSubset(A,B)){
					if(simplex->weight < y->weight)
						simplex->weight = y->weight;
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
int di=0;
for( auto x : this->simplexList)
	std::cout<<"Count of "<<di++<<"-simplex ::"<<x.size()<<"\n";
	return;
}




template<typename nodeType>
void alphaComplex<nodeType>::buildBetaComplexFilteration(std::vector<std::vector<unsigned>> dsimplexmesh, int npts, std::vector<std::vector<double>> inputData,kdTree tree){	
	this->ut.writeLog(this->simplexType,"No build alpha complex function defined");
	return;
}

template<typename nodeType>
void alphaComplex<nodeType>::buildFilteration(std::vector<std::vector<unsigned>> dsimplexmesh, int npts, std::vector<std::vector<double>> inputData,double beta,kdTree){	
	this->ut.writeLog(this->simplexType,"No build beta complex function defined");
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

