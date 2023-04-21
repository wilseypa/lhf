/**
 * @file utils.hpp
 *
 * @brief Contains utility functions for the LFH system (https://github.com/wilseypa/LFH).
 */

#include <string>
#include <cmath>
#include <algorithm>
#include <utility>
#include <numeric>
#include <iostream>
#include <fstream>
#include "utils.hpp"
#include <time.h>

/**
 * @brief Contains utility functions for the LFH system (https://github.com/wilseypa/LFH).
 */

 /**
     * @brief Extracts boundary points from a vector of simplex nodes.
     *
     * @tparam T The type of the simplex nodes.
     * @param nodes The vector of simplex nodes.
     * @return A set of indices corresponding to the boundary points.
     */

template std::set<unsigned int, std::less<unsigned int>, std::allocator<unsigned int> > utils::extractBoundaryPoints<simplexNode>(std::vector<std::shared_ptr<simplexNode>, std::allocator<std::shared_ptr<simplexNode> > >);
template std::set<unsigned int, std::less<unsigned int>, std::allocator<unsigned int> > utils::extractBoundaryPoints<alphaNode>(std::vector<std::shared_ptr<alphaNode>, std::allocator<std::shared_ptr<alphaNode> > >);
template std::set<unsigned int, std::less<unsigned int>, std::allocator<unsigned int> > utils::extractBoundaryPoints<witnessNode>(std::vector<std::shared_ptr<witnessNode>, std::allocator<std::shared_ptr<witnessNode> > >);

/**
 * @brief Initializes a Union-Find data structure.
 *
 * @param n The number of elements in the Union-Find data structure.
 */

unionFind::unionFind(int n) : rank(n, 0), parent(n, 0) {
	for(int i=0; i<n; i++) parent[i]=i;
}

/**
  @file utils.hpp
 
  @brief Contains utility functions for the LFH system (https://github.com/wilseypa/LFH).
 */

/**
  @brief Contains utility functions for the LFH system (https://github.com/wilseypa/LFH).
 */

/**
      @brief Extracts boundary points from a vector of simplex nodes.
     
      @tparam T The type of the simplex nodes.
      @param nodes The vector of simplex nodes.
      @return A set of indices corresponding to the boundary points.
     */

	/**
  @brief Finds the root of the component that contains the given element.
 
  @param i The index of the element.
  @return The index of the root of the component.
 */

int unionFind::find(int i){
	if(i == parent[i]) return i; //Found name of the component
	parent[i] = find(parent[i]); //Path Compression
	return parent[i];
}

/**
  @brief Joins two components together.
 
  @param x The index of an element in one component.
  @param y The index of an element in another component.
  @return True if the two components were joined, false otherwise.
 */

bool unionFind::join(int x, int y){ //Union by rank
	x = find(x);
	y = find(y);
	if(x == y) return false;
	if(rank[x] == rank[y]){
		rank[y]++;
		parent[x] = y;
	} else if(rank[x] < rank[y]){
		parent[x] = y;
	} else{
		parent[y] = x;
	}
	return true;
}

/**
  @brief Contains utility functions for the LFH system (https://github.com/wilseypa/LFH).
 */

/**
      @brief Determines whether A is a subset of B.
     
      @param A The vector to be tested for being a subset.
      @param B The vector that A is being tested against.
      @return True if A is a subset of B, false otherwise.
     */

bool utils:: isSubset(std::vector<unsigned> A,std::vector<unsigned> B){
	std::sort(A.begin(),A.end());
	std::sort(B.begin(),B.end());
	return std::includes(A.begin(),A.end(),B.begin(),B.end());
}

  /**
      @brief Constructs a new instance of the utils class.
     */

// utils constructor, currently no needed information for the class constructor
utils::utils(){}

/**
      @brief Constructs a new instance of the utils class.
     
      @param _debug The debug level.
      @param _outputFile The output file.
	 */

utils::utils(std::string _debug, std::string _outputFile){
	debug = _debug;
	outputFile = _outputFile;
}

/**
      @brief Calculates the average of a vector of doubles.
     
      @param v The vector of doubles.
      @return The average of the vector of doubles.
     */

/**
  @brief Determines whether A is a subset of B.
 
  @param A The vector to be tested for being a subset.
  @param B The vector that A is being tested against.
  @return True if A is a subset of B, false otherwise.
 */

/**
  @brief Constructs a new instance of the utils class.
 */

/**
  @brief Constructs a new instance of the utils class.
 
  @param _debug The debug level.
  @param _outputFile The output file.
 */

/**
  @brief Calculates the average of a vector of doubles.
 
  @param v The vector of doubles.
  @return The average of the vector of doubles.
 */

double  utils :: getAverage(std::vector<double> &v) {
    if (v.empty()) {
        return 0;
    }
    return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
}

std::vector<std::vector<double>> utils ::  transposeMeanAdjusted(std::vector<std::vector<double>> input){
    std::vector<std::vector<double>> inputtranspose;
    std::vector<std::vector<double>> inputstd;
  
    for(int i=0;i<input[0].size();i++){
		std::vector<double> row;
	  for(int j=0;j<input.size();j++){
		  row.push_back(input[j][i]);
	}
	inputtranspose.push_back(row);
	}
	
    for(int i =0;i<inputtranspose.size();i++){
		std::vector<double> inter;
		double averagex = getAverage(inputtranspose[i]);
		double stdval = 0;
		for(int j =0;j<inputtranspose[i].size();j++){
			inter.push_back((inputtranspose[i][j] -averagex));
		}
		inputstd.push_back(inter);
		}
	return inputstd;
}

Eigen::MatrixXd utils ::  covariance(std::vector<std::vector<double>> input){
	Eigen::MatrixXd result(input.size(),input.size());
	int i=0,j=0;
	for(auto x : input){
		for(auto y : input){
			double value = 0;
			double averagex = getAverage(x);
			double averagey = getAverage(y);
			for(int i =0;i<x.size();i++){
				value += (x[i] -averagex)*(y[i] -averagey);
			}
			result(i,j) = value/(x.size()-1);
			j++;
		}
		i++;
		j=0;
	}
	return result;
}
	

double utils :: determinantOfMatrix(std::vector<std::vector<double>> mat, int n)
{
  double  det = 1;
  int index;
  for (unsigned i = 0; i < n; i++){
        index = i;
				while (mat[index][i] == 0 && index < n)
				       index++;
	      if (index == n)
				       continue;
				if (index != i){
				  for (int j = 0; j < n; j++){
						  double temp12 = mat[index][j];
				      mat[index][j] = mat[i][j];
				      mat[i][j] = temp12;
				  }
				  det = det * pow(-1, index - i);
			  }
        double rectemp = mat[i][i];
        for (unsigned j = i; j < n; j++)
                mat[i][j] /= rectemp;
        for (unsigned j = i + 1; j < n; j++){
					if(mat[j][i] != 0){
					      double rectemp2 = mat[j][i];
                for (unsigned t = i;  t< n; t++)
									             mat[i][t] *= rectemp2;
                for (unsigned t = i;  t< n; t++)
                               mat[j][t] -= mat[i][t];
                for (unsigned k = 0; k < n; k++)
                           		 mat[i][k] /= rectemp2;
						}
        }
        for (unsigned k = 0; k < n; k++)
              mat[i][k] *= rectemp;
      }
  for (unsigned i = 0; i < n; i++)
	        det = det * mat[i][i];

	return (det);
}

std::vector<std::vector<double>> utils ::  transpose(std::vector<std::vector<double>> input){
    std::vector<std::vector<double>> inputtranspose;
    
    for(int i=0;i<input[0].size();i++){
		std::vector<double> row;
	  for(int j=0;j<input.size();j++){
		  row.push_back(input[j][i]);
	}
	inputtranspose.push_back(row);
	}
	return inputtranspose;
}


std::vector<std::vector<double>> utils :: matrixMultiplication(std::vector<std::vector<double>> matA, std::vector<std::vector<double>> matB){
			int n1 = matA.size();
			int m1 = matA[0].size();
			int n2 = matB.size();
			int m2 = matB[0].size();
			std::vector<std::vector<double>>  mat(n1,std::vector<double>(m2,0));

			if(m1!=n2)
					return mat;

	    for (int i = 0; i < n1; i++){
	        for (int j = 0; j < m2; j++){
	            for (int x = 0; x < m1; x++){
	                mat[i][j] += matA[i][x]* matB[x][j];
	            }
	        }
	    }
	return mat;
}

std::vector<std::vector<double>> utils :: inverseOfMatrix(std::vector<std::vector<double>> mat, int n){
    int index;
    std::vector<std::vector<double>> matinv(n,std::vector<double> (n,0));

    for(int i =0 ;i<n;i++)
      matinv[i][i] = 1;

    for (unsigned i = 0; i < n; i++)
    {
        index = i;
        while (mat[index][i] == 0 && index < n)
      			index++;
        if (index == n)
            continue;
        if (index != i){
            for (int j = 0; j < n; j++){
                double temp12 = mat[index][j];
                  mat[index][j] = mat[i][j];
                  mat[i][j] = temp12;
                double temp121 = matinv[index][j];
                    matinv[index][j] = matinv[i][j];
                    matinv[i][j] = temp121;
            }
				}
        double rectemp = mat[i][i];
        if(mat[i][i]!=1){
          for (unsigned j = 0; j < n; j++){
                mat[i][j] /= rectemp;
                matinv[i][j] /= rectemp;
              }
        }
        for (unsigned j = 0; j < n; j++)
        {
           if(mat[j][i] != 0 && j!=i){
           			double rectemp2 = mat[j][i];
           			for (unsigned t = 0;  t< n; t++){
             				mat[i][t] *= rectemp2;
             				matinv[i][t] *= rectemp2;
            		}
            		for (unsigned t = 0;  t< n; t++){
              			mat[j][t] -= mat[i][t];
              			matinv[j][t] -= matinv[i][t];
             		}
               for (unsigned k = 0; k < n; k++){
                   mat[i][k] /= rectemp2;
                   matinv[i][k] /= rectemp2;
								}
					}
			}
	}
  return matinv;
}

std::vector<std::vector<double>> utils :: betaCentersCalculation(std::vector<double> hpcoff, double beta, double circumRadius,std::vector<double> circumCenter){
	double distance = sqrt(pow((beta*circumRadius),2) - pow(circumRadius,2));
	double d1 , d2;   // Parallel Plane coefficient
	double sqrtofsquaredsum =0,squaredsum=0;
	double dotproduct = 0;
    int i=0;
	for(auto x: hpcoff){
		squaredsum += x*x;
	    dotproduct += x*circumCenter[i];
	    i++;
	}
	
	sqrtofsquaredsum = sqrt(squaredsum);
	 
	d1 = -dotproduct + distance*sqrtofsquaredsum;
	d2 = -dotproduct - distance*sqrtofsquaredsum; 
	
	double t1 , t2;
	
	t1 = (-dotproduct -d1)/squaredsum;
	
	t2 = (-dotproduct -d2)/squaredsum;
	
	std::vector<std::vector<double>> centers;
	std::vector<double> center1;
	std::vector<double> center2;
	i=0;
	for(auto x: hpcoff){
	    center1.push_back(x*t1 + circumCenter[i]);
	    center2.push_back(x*t2 + circumCenter[i]);
	    i++;
	}
	centers.push_back(center1);
	centers.push_back(center2);
	return centers;
}
std::pair<std::vector<std::vector<double>>,std::vector<std::vector<double>>> utils:: computePCA(std::vector<std::vector<double>> input, int targetDim){
        std::vector<std::vector<double>> transp = transposeMeanAdjusted(input);
        Eigen::MatrixXd covar = covariance(transp);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
        es.compute(covar);		
		std::vector<double> eigenval;
	    std::vector<std::vector<double>> eigenvec;
		for(int k=input[0].size()-1;k>=0;k--){
			std::vector<double> evec;	
			for(int p=0;p<input[0].size();p++)
				evec.push_back(es.eigenvectors().col(k)[p]);
			eigenvec.push_back(evec);
		}	

        for(int p=input[0].size()-1;p>=0;p--)
			eigenval.push_back(es.eigenvalues()[p]);
	    std::vector<std::vector<double>> eigenVecSelect;

		for(int k=0;k<targetDim;k++)
			eigenVecSelect.push_back(eigenvec[k]);
		
		return std::make_pair(transpose(matrixMultiplication(eigenVecSelect,transp)),eigenVecSelect);

		
}
std::vector<std::vector<double>> utils:: computePCAInverse(std::vector<std::vector<double>> input,std::vector<std::vector<double>> FinalOutput, std::vector<std::vector<double>> eigenvectors ){
	    
	    auto result = transpose(matrixMultiplication(transpose(eigenvectors),transpose(FinalOutput)));
	    
	    std::vector<std::vector<double>> transposeinput = transpose(input);
	    int i =0;
	    for(auto x:transposeinput){
		    double average = getAverage(x);
		    for(int y=0;y<result.size();y++){
				result[y][i] +=average;
		     }
		    i++;
		}
        return result;

}
std::pair<std::vector<double>,std::vector<std::vector<double>>> utils :: nullSpaceOfMatrix(std::set<unsigned> simplex,std::vector<std::vector<double>> inputData,std::vector<double> & cc, double radius,bool lowerdimension){
    int index;
    srand(time(NULL));
    int n = simplex.size();
    std::vector<double> matns(n,1);
    std::vector<std::vector<double>> mat;
    
    for(auto x : simplex)
	    mat.push_back(inputData[x]);
	mat.push_back(cc);
	std::pair<std::vector<std::vector<double>>,std::vector<std::vector<double>>> outputPCA;
	if(lowerdimension==true){    
		outputPCA = computePCA(mat,simplex.size());
		mat = outputPCA.first;
		cc = mat[mat.size()-1];
		mat.pop_back();
    }
   	
    for (unsigned i = 0; i < n; i++)
    {
        index = i;
        while (mat[index][i] == 0 && index < n)
      			index++;
        if (index == n)
            continue;
        if (index != i){
            for (int j = 0; j < n; j++){
                double temp12 = mat[index][j];
                  mat[index][j] = mat[i][j];
                  mat[i][j] = temp12;
			  }
                double temp121 = matns[index];
                    matns[index] = matns[i];
                    matns[i] = temp121;
            
		}
        double rectemp = mat[i][i];
        if(mat[i][i]!=1){
          for (unsigned j = 0; j < n; j++){
                mat[i][j] /= rectemp;
			}
                matns[i] /= rectemp;
        }
     
        for (unsigned j = 0; j < n; j++)
        {
           if(mat[j][i] != 0 && j!=i){
           			double rectemp2 = mat[j][i];
           			for (unsigned t = 0;  t< n; t++){
             				mat[i][t] *= rectemp2;
             		}
             				matns[i] *= rectemp2;
            		for (unsigned t = 0;  t< n; t++){
              			mat[j][t] -= mat[i][t];
              		}
              			matns[j] -= matns[i];
               for (unsigned k = 0; k < n; k++){
                   mat[i][k] /= rectemp2;
			   }
                   matns[i] /= rectemp2;
			}
		}
 }

  return std::make_pair(matns,outputPCA.second);
}
std::vector<double> utils :: circumCenter(std::set<unsigned> simplex,std::vector<std::vector<double>> inputData){
// Soluiton  = inv(matA) * matC
   std::vector<std::vector<double>>  matA(simplex.size());
	 std::vector<std::vector<double>>  invmatA;
	 std::vector<std::vector<double>>  matC(simplex.size());
	 std::vector<std::vector<double>> rawCircumCenter;
	 std::vector<double> circumCenter;
	 std::set<unsigned> simplexcopy = simplex;

	 auto it = simplex.end();
	 it--;
	 int ii =0;
   unsigned Sn = *(it);
	 simplex.erase(Sn);
	 for(auto i : simplex){
	 		for(auto j : simplex){
				  std::vector<double> d1,d2;
					double dotProduct=0;
					std::transform(inputData[i].begin(), inputData[i].end(), inputData[Sn].begin(), std::back_inserter(d1),[](double e1,double e2){return (e1-e2);});
					std::transform(inputData[j].begin(), inputData[j].end(), inputData[Sn].begin(), std::back_inserter(d2),[](double e1,double e2){return (e1-e2);});
					for (int k = 0; k < inputData[0].size(); k++){
              dotProduct = dotProduct + d1[k] * d2[k];
						}
					matA[ii].push_back(dotProduct);
				 if(i==j)
				 	matC[ii].push_back(dotProduct/2);
			}
			matA[ii].push_back(0);
			ii++;
		}
		for(int i =0;i<simplex.size()+1;i++)
			matA[simplex.size()].push_back(1);
		matC[simplex.size()].push_back(1);
		invmatA = inverseOfMatrix(matA,matA[0].size());
		rawCircumCenter = matrixMultiplication(invmatA,matC);
		for(int i = 0; i < inputData[0].size(); i++){
       double coordinate = 0;
			 std::set<unsigned> ::iterator index = simplexcopy.begin();
			for(int j =0;j<rawCircumCenter.size();j++){
				coordinate += rawCircumCenter[j][0]*inputData[(*index)][i];
				index++;
			}
			circumCenter.push_back(coordinate);
		}



		return circumCenter;

}
double utils :: circumRadius(std::set<unsigned> simplex,std::vector<std::vector<double>>* distMatrix){
    std::vector<std::vector<double>>  matA(simplex.size());
		std::vector<std::vector<double>>  matACap(simplex.size()+1);
		int ii=0;
	  for(auto i : simplex){
			matACap[ii+1].push_back(1);
			for(auto j :simplex){
				if((*distMatrix)[i][j]!=0){
		   	matA[ii].push_back(pow(((*distMatrix)[i][j]),2));
				matACap[ii+1].push_back(pow(((*distMatrix)[i][j]),2));
			}
			else{
				matA[ii].push_back(pow(((*distMatrix)[j][i]),2));
		  	matACap[ii+1].push_back(pow(((*distMatrix)[j][i]),2));
			}
	  }
		ii++;
	}
	matACap[0].push_back(0);
	for(auto i : simplex)
    matACap[0].push_back(1);

	return -(determinantOfMatrix(matA,simplex.size())/(2*determinantOfMatrix(matACap,simplex.size()+1)));
}
double utils :: simplexVolume(std::set<unsigned> simplex,std::vector<std::vector<double>>* distMatrix,int dd){
		std::vector<std::vector<double>>  matACap(simplex.size()+1);
		int ii=0;
	  for(auto i : simplex){
			matACap[ii+1].push_back(1);
			for(auto j :simplex){
				if((*distMatrix)[i][j]!=0){
				matACap[ii+1].push_back(pow(((*distMatrix)[i][j]),2));
			}
			else{
		  	matACap[ii+1].push_back(pow(((*distMatrix)[j][i]),2));
			}
	  }
		ii++;
	}
	matACap[0].push_back(0);
	for(auto i : simplex)
		matACap[0].push_back(1);
    if(dd%2==0)
		return ((-1)*determinantOfMatrix(matACap,simplex.size()+1))/(pow(2,dd)*(pow(tgamma(dd+1),2)));
	else
		return (determinantOfMatrix(matACap,simplex.size()+1))/(pow(2,dd)*(pow(tgamma(dd+1),2)));
}
double utils :: simplexVolume(std::vector<std::vector<double>> spoints){
		std::vector<std::vector<double>>  matACap(spoints.size()+1);
		int ii=0;
	  for(auto i : spoints){
			matACap[ii+1].push_back(1);
			for(auto j :spoints){
				matACap[ii+1].push_back(pow(vectors_distance(i,j),2));
			}
		ii++;
	}
	matACap[0].push_back(0);
	for(int i =0; i< spoints.size(); i++)
		matACap[0].push_back(1);
    if(spoints.size()%2==0)
		return ((-1)*determinantOfMatrix(matACap,spoints.size()+1))/(pow(2,spoints[0].size())*(pow(tgamma(spoints[0].size()+1),2)));
	else
		return (determinantOfMatrix(matACap,spoints.size()+1))/(pow(2,spoints[0].size())*(pow(tgamma(spoints[0].size()+1),2)));
}

double utils::computeMaxRadius(int k, std::vector<std::vector<double>> &centroids, std::vector<std::vector<double>> &originalData, std::vector<unsigned> &labels){
	double maxRadius = 0;
	double curRadius = 0;

	//Iterate through each point
	for(unsigned i = 0; i < originalData.size(); i++){
		//Check the distance of this point to it's centroid
		curRadius = vectors_distance(originalData[i], centroids[labels[i]]);

		if(curRadius > maxRadius)
			maxRadius = curRadius;
	}

	return maxRadius;
}

double utils::computeAvgRadius(int k, std::vector<std::vector<double>> &centroids, std::vector<std::vector<double>> &originalData, std::vector<unsigned> &labels){
	double totalRadius = 0;

	//Iterate through each point
	for(unsigned i = 0; i < originalData.size(); i++){

		//Check the distance of this point to it's centroid
		totalRadius += vectors_distance(originalData[i], centroids[labels[i]]);
	}

	return totalRadius / originalData.size();
}

std::vector<std::vector<bool>> utils :: betaNeighbors(std::vector<std::vector<double>>& inData,double beta,std::string betaMode){
	std::vector<std::vector<bool>> incidenceMatrix(inData.size(),std::vector<bool>(inData.size(),0));
	kdTree tree(inData, inData.size()); //KDTree for efficient nearest neighbor search
    double betafactor;
    std::vector<std::vector<double>> rotationMatrix(inData[0].size(),std::vector<double>(inData[0].size(),0));

//[
//cos(90) -sin(90)  ................... 0
//sin(90) cos(90)  .....................0
//0		0	0 ......0.......0
//......................1......	0	0
//0		0	0	1	0
//0		0	0	0    	1
//]


    for(unsigned i =0;i<inData[0].size();i++){
		for(unsigned j = 0;j<inData[0].size();j++){
			if(i==1 && j==0)
				rotationMatrix[i][j] = 1;
			else if(i==0 && j==1)
			    rotationMatrix[i][j] = -1;
			else if(i==0 && j==0)
			     rotationMatrix[i][j] =0;
			else if(i==1 && j==1)
			     rotationMatrix[i][j] =0;
			else if(i==j)
				 rotationMatrix[i][j] = 1;
			else
				 rotationMatrix[i][j] = 0;
		}
	}

	if(beta < 0){
		std::cout<<"Invalid Beta Value :: Value set to 1"<<std::endl;
		beta = 1;  // Changing to Default 
	}	
	for(unsigned i = 0; i < inData.size(); i++){
		for(unsigned j = i+1; j < inData.size(); j++){
		/*	std::cout<<" First Point \n";
			for(auto x:inData[i])
				std::cout<<x<<" ";

			std::cout<<" Second Point \n";
			for(auto x:inData[j])
				std::cout<<x<<" ";
		*/	double radius;	
			std::vector<double> temp;
			std::vector<double> center1;
			std::vector<double> center2;
			std::vector<double> temp1;
			std::vector<double> tempfinal;
			std::vector<std::vector<double>> colvector;
			if(beta <=1){   // Common Definition for lune and circle based beta-Skeleton
		              radius = vectors_distance(inData[i],inData[j])/(2*beta);
		              std::transform(inData[i].begin(), inData[i].end(), inData[j].begin(), std::back_inserter(temp),[](double e1, double e2) {return ((e1+e2)/2);});
		              betafactor = sqrt((1-pow(beta,2)))/(2*beta);
		              std::transform(inData[i].begin(), inData[i].end(), inData[j].begin(), std::back_inserter(temp1),[](double e1, double e2) {return ((e2-e1));});
                      for(auto x : temp1){
						  std::vector<double> singlevector;
						  singlevector.push_back(x);
						  colvector.push_back(singlevector);
					  }
					  std::vector<std::vector<double>> matmul = matrixMultiplication(rotationMatrix,colvector);
					  for(auto x: matmul)
						for(auto y : x)
							tempfinal.push_back(y*betafactor);	
					   std::transform(temp.begin(), temp.end(), tempfinal.begin(), std::back_inserter(center1),[](double e1, double e2) {return ((e1+e2));});
                      			   std::transform(temp.begin(), temp.end(), tempfinal.begin(), std::back_inserter(center2),[](double e1, double e2) {return ((e1-e2));});
                       
                       			   
			/*	     	   std::cout<<" First Center and radius:: "<<radius<<"\n";
				  	   for(auto x:center1)
						std::cout<<x<<" ";
				     	   std::cout<<" Second Center \n";
				  	   for(auto x:center2)
						std::cout<<x<<" ";
			*/		   std::vector<size_t> neighbors1 = tree.neighborhoodIndices(center1, radius); //All neighbors in radius-ball
					   std::vector<unsigned> toremove;
					   toremove.push_back(i);
					   toremove.push_back(j); 
					   std::sort(toremove.begin(),toremove.end());
					   std::sort(neighbors1.begin(), neighbors1.end());
					   std::vector<int> neg1;
					   std::set_difference(neighbors1.begin(),neighbors1.end(),toremove.begin(),toremove.end(), std::back_inserter(neg1));
		    
    		
					   std::vector<size_t> neighbors2 = tree.neighborhoodIndices(center2, radius); //All neighbors in radius-ball
					   std::sort(neighbors2.begin(), neighbors2.end());
					   std::vector<int> neg2;
					   std::set_difference(neighbors2.begin(),neighbors2.end(),toremove.begin(),toremove.end(), std::back_inserter(neg2));
		    
    		
    		
					   std::vector<size_t> neighbors;
          
					   std::vector<size_t> v(std::min(neg1.size(),neg2.size()));
					   std::vector<size_t>::iterator it;
					   std::sort (neg1.begin(),neg1.end());
					   std::sort (neg2.begin(),neg2.end()); 
					   it=std::set_intersection (neg1.begin(), neg1.end(), neg2.begin(), neg2.end(), v.begin());
					   v.resize(it-v.begin()); 
					   neighbors = v;				
		
					   if(neighbors.size() == 0)
						 incidenceMatrix[i][j] = true;
					   else
						 incidenceMatrix[i][j] = false;
				
			}
			else if(betaMode == "lune"){ // Lune Based Beta Skeleton for beta > 1
			          radius = (vectors_distance(inData[i],inData[j])*beta)/2;
                      double bf = beta/2;
                      double bf1 = 1-beta/2;

                      std::transform(inData[i].begin(), inData[i].end(), inData[j].begin(), std::back_inserter(center1),[bf,bf1](double e1, double e2) {return ((e1*bf+e2*bf1));});
					  std::transform(inData[j].begin(), inData[j].end(), inData[i].begin(), std::back_inserter(center2),[bf,bf1](double e1, double e2) {return ((e1*bf+e2*bf1));});
					   
					   std::vector<size_t> neighbors1 = tree.neighborhoodIndices(center1, radius); //All neighbors in radius-ball
					   std::vector<unsigned> toremove;
					   toremove.push_back(i);
					   toremove.push_back(j); 
					   std::sort(toremove.begin(),toremove.end());
					   std::sort(neighbors1.begin(), neighbors1.end());
					   std::vector<int> neg1;
					   std::set_difference(neighbors1.begin(),neighbors1.end(),toremove.begin(),toremove.end(), std::back_inserter(neg1));
		    
					   std::vector<size_t> neighbors2 = tree.neighborhoodIndices(center2, radius); //All neighbors in radius-ball
					    std::sort(neighbors2.begin(), neighbors2.end());
					   std::vector<int> neg2;
					   std::set_difference(neighbors2.begin(),neighbors2.end(),toremove.begin(),toremove.end(), std::back_inserter(neg2));
		    
		    
					   std::vector<size_t> neighbors;
          
					   std::vector<size_t> v(std::min(neg1.size(),neg2.size()));
					   std::vector<size_t>::iterator it;
					   std::sort (neg1.begin(),neg1.end());
					   std::sort (neg2.begin(),neg2.end()); 
					   it=std::set_intersection (neg1.begin(), neg1.end(), neg2.begin(), neg2.end(), v.begin());
					   v.resize(it-v.begin()); 
					   neighbors = v;				
		
					   if(neighbors.size() == 0)
						 incidenceMatrix[i][j] = true;
					   else
						 incidenceMatrix[i][j] = false;
					  	 
			}
			else if(betaMode == "circle"){			// Circle Based Beta Skeleton for beta >1
					  radius = (vectors_distance(inData[i],inData[j])*beta)/(2);
		              std::transform(inData[i].begin(), inData[i].end(), inData[j].begin(), std::back_inserter(temp),[](double e1, double e2) {return ((e1+e2)/2);});
		              betafactor = sqrt((pow(beta,2)-1))/2;
		              std::transform(inData[i].begin(), inData[i].end(), inData[j].begin(), std::back_inserter(temp1),[](double e1, double e2) {return ((e2-e1));});
		              for(auto x : temp1){
						  std::vector<double> singlevector;
						  singlevector.push_back(x);
						  colvector.push_back(singlevector);
					  }
					  std::vector<std::vector<double>> matmul = matrixMultiplication(rotationMatrix,colvector);
					  for(auto x: matmul)
						for(auto y : x)
							tempfinal.push_back(y*betafactor);	
                      std::transform(temp.begin(), temp.end(), tempfinal.begin(), std::back_inserter(center1),[](double e1, double e2) {return ((e1+e2));});
                      std::transform(temp.begin(), temp.end(), tempfinal.begin(), std::back_inserter(center2),[](double e1, double e2) {return ((e1-e2));});
                       
					  std::vector<size_t> neighbors1 = tree.neighborhoodIndices(center1, radius); //All neighbors in radius-ball
					   std::vector<unsigned> toremove;
					   toremove.push_back(i);
					   toremove.push_back(j); 
					   std::sort(toremove.begin(),toremove.end());
					   std::sort(neighbors1.begin(), neighbors1.end());
					   std::vector<int> neg1;
					   std::set_difference(neighbors1.begin(),neighbors1.end(),toremove.begin(),toremove.end(), std::back_inserter(neg1));
		    
					  std::vector<size_t> neighbors2 = tree.neighborhoodIndices(center2, radius); //All neighbors in radius-ball
					   std::sort(neighbors2.begin(), neighbors2.end());
					   std::vector<int> neg2;
					   std::set_difference(neighbors2.begin(),neighbors2.end(),toremove.begin(),toremove.end(), std::back_inserter(neg2));
		    
					  std::vector<size_t> neighbors;
	            
					  std::vector<size_t> v(neighbors1.size() + neighbors2.size());
					  std::vector<size_t>::iterator it;	
					  std::sort (neg1.begin(),neg1.end());
					  std::sort (neg2.begin(),neg2.end()); 
					  it=std::set_union (neg1.begin(), neg1.end(), neg2.begin(), neg2.end(), v.begin());
					  v.resize(it-v.begin()); 
					  neighbors = v;
			          if(neighbors.size() == 0)
						 incidenceMatrix[i][j] = true;
					   else
						 incidenceMatrix[i][j] = false;
						 
       			}
		}
	}
	return incidenceMatrix;
}





std::pair<std::vector<std::vector<unsigned>>, std::vector<std::vector<std::vector<double>>>> utils::separatePartitions(int k, std::vector<std::vector<double>> originalData, std::vector<unsigned> labels){
	std::vector<std::vector<double>> a;
	std::vector<unsigned> b;
	std::vector<std::vector<std::vector<double>>> res(k, a);
	std::vector<std::vector<unsigned>> labres(k, b);

	for(unsigned i = 0; i < labels.size(); i++){
		res[labels[i]].push_back(originalData[i]);
		labres[labels[i]].push_back(i);
	}

	return std::make_pair(labres, res);
}

std::pair<std::vector<std::vector<unsigned>>, std::vector<std::vector<std::vector<double>>>> utils::separatePartitions(double rad, std::vector<std::vector<double>> centroids, std::vector<std::vector<double>> originalData, std::vector<unsigned> labels){
	std::vector<std::vector<double>> a;
	std::vector<unsigned> b;

	//Store the results to return
	std::vector<std::vector<std::vector<double>>> res(centroids.size(), a);
	std::vector<std::vector<unsigned>> labres(centroids.size(), b);

	//Iterate through each label
	for(unsigned i = 0; i < labels.size(); i++){

		//Check for this point belonging to each centroid
		for(unsigned j = 0; j < centroids.size(); j++){
			if(labels[i] == j){
				//If this is a labeled constituent put to front of array
				res[j].insert(res[j].begin(), originalData[i]);
				labres[j].insert(labres[j].begin(), i);

			}else{

				double curRad = vectors_distance(originalData[i], centroids[j]);

				//If this distance is less than our radius cutoff push to back
				if(curRad < rad){
					res[j].push_back(originalData[i]);
					labres[j].push_back(i);
				}
			}
		}
	}

	return std::make_pair(labres, res);
}

std::vector<std::vector<std::vector<double>>> utils::separateBoundaryPartitions(std::vector<std::set<unsigned>> boundaryLists, std::vector<std::vector<double>> originalData, std::vector<unsigned> labels){
	std::vector<std::vector<double>> a;
	std::vector<std::vector<std::vector<double>>> res(boundaryLists.size(), a);

	for(unsigned i = 0; i < originalData.size(); i++){

		for(unsigned j = 0; j < boundaryLists.size(); j++){

			if(boundaryLists[j].find(labels[i]) != boundaryLists[j].end())
				res[j].push_back(originalData[i]);
		}
	}

	return res;
}

// void utils::extractBoundaryPoints(std::vector<bettiBoundaryTableEntry>& bettiTable){
// 	for(auto& bet: bettiTable){
// 		std::set<unsigned> bound;
// 		for(auto simplex : bet.boundary) bound.insert(simplex->simplex.begin(), simplex->simplex.end());
// 		bet.boundaryPoints = bound;
// 	}
// }

template <typename T>
std::set<unsigned> utils::extractBoundaryPoints(std::vector<std::shared_ptr<T>> boundary){
	std::set<unsigned> boundaryPoints;
	for(auto simplex : boundary) boundaryPoints.insert(simplex->simplex.begin(), simplex->simplex.end());
	return boundaryPoints;
}

template <typename T>
std::set<unsigned> utils::extractBoundaryPoints(std::vector<T*> boundary){
	std::set<unsigned> boundaryPoints;
	for(auto simplex : boundary) boundaryPoints.insert(simplex->simplex.begin(), simplex->simplex.end());
	return boundaryPoints;
}

std::vector<bettiBoundaryTableEntry> utils::mapPartitionIndexing(std::vector<unsigned> partitionedLabels, std::vector<bettiBoundaryTableEntry> bettiTable){
	for(auto& bet : bettiTable){
		std::set<unsigned> convBound;

		for(auto ind : bet.boundaryPoints){
			convBound.insert(partitionedLabels[ind]);
		}

		bet.boundaryPoints = convBound;
	}
	return bettiTable;
}

void utils::print2DVector(const std::vector<std::vector<unsigned>>& a){
	for(unsigned i = 0; i < a.size(); i++){
		for(unsigned j = 0; j < a[i].size(); j++){
			std::cout << a[i][j] << '\t';
		}
		std::cout << std::endl;
	}
	return;
}

void utils::print1DSet(const std::pair<std::set<unsigned>, double>& a){
		std::cout << "Test\t";

		for(auto iter = a.first.begin(); iter!= a.first.end(); iter++){
			std::cout << *iter << ",";
		}
		std::cout << "\t";
}

void utils::print1DVector(const std::vector<double>& a){
	for(unsigned i = 0; i < a.size(); i++){
			std::cout << a[i] << ",";
	}
	std::cout << "\n";
	return;
}

void utils::print1DVector(const std::vector<unsigned>& a){
	for(unsigned i = 0; i < a.size(); i++){
			std::cout << a[i] << ",";
	}
	std::cout << "\n";
	return;
}

void utils::print1DVector(const std::set<unsigned>& a){
	for(auto z : a){
			std::cout << z << ",";
	}
	std::cout << "\n";
	return;
}

double utils::vectors_distance(const double& a, const double& b){
		return pow((a-b),2);
}

std::set<unsigned> utils::setXOR(std::set<unsigned>& setA, std::set<unsigned>& setB){
	std::set<unsigned> ret;

	//if(setA.size() == 0)
	//	ret = setB;
	//else if (setB.size() == 0)
	//	ret = setA;
	//else
	set_symmetric_difference(setA.begin(), setA.end(), setB.begin(), setB.end(), std::inserter(ret, ret.begin()));

	return ret;
}

double utils::vectors_distance(const std::vector<double>& a, const std::vector<double>& b){
	std::vector<double> temp;

	if(b.size() == 0)
		return 0;

	std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(temp),[](double e1, double e2) {
		return pow((e1-e2),2);
	});

	return sqrt(std::accumulate(temp.begin(), temp.end(), 0.0));
}

std::vector<unsigned> utils::setIntersect(std::vector<unsigned> v1, std::vector<unsigned> v2, bool isSorted){
	std::vector<unsigned> ret;

	if(v1 == v2)
		return v1;

	if(!isSorted){
		sort(v1.begin(), v1.end());
		sort(v2.begin(), v2.end());
	}

	set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(ret));

	/*for(auto iter = v1.begin(); iter!= v1.end(); iter++){
		std::cout << *iter << ",";
	}
	std::cout << "\t";
	for(auto iter = v2.begin(); iter!= v2.end(); iter++){
		std::cout << *iter << ",";
	}
	std::cout << "\t";
	for(auto iter = ret.begin(); iter!= ret.end(); iter++){
		std::cout << *iter << ",";
	}
	std::cout << std::endl;*/

	return ret;

}

std::set<unsigned> utils::setIntersect(std::set<unsigned> v1, std::set<unsigned> v2, bool isSorted){
	std::set<unsigned> ret;

	if(v1 == v2)
		return v1;

	//if(!isSorted){
	//	sort(v1.begin(), v1.end());
	//	sort(v2.begin(), v2.end());
	//}

	set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), std::inserter(ret, ret.begin()));

	/*for(auto iter = v1.begin(); iter!= v1.end(); iter++){
		std::cout << *iter << ",";
	}
	std::cout << "\t";
	for(auto iter = v2.begin(); iter!= v2.end(); iter++){
		std::cout << *iter << ",";
	}
	std::cout << "\t";
	for(auto iter = ret.begin(); iter!= ret.end(); iter++){
		std::cout << *iter << ",";
	}
	std::cout << std::endl;*/

	return ret;

}



// Find the intersect of two vectors
std::pair<std::vector<unsigned>, std::vector<unsigned>> utils::intersect(std::vector<unsigned> v1, std::vector<unsigned> v2, bool isSorted){
	std::pair<std::vector<unsigned>, std::vector<unsigned>> ret;
	std::pair<std::vector<unsigned>, std::vector<unsigned>> retTemp;

	if(v1 == v2)
		return ret;

	if(!isSorted){
		sort(v1.begin(), v1.end());
		sort(v2.begin(), v2.end());
	}

	set_union(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(retTemp.second));
	set_symmetric_difference(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(retTemp.first));

	/*for(auto iter = v1.begin(); iter!= v1.end(); iter++){
		std::cout << *iter << ",";
	}
	std::cout << "\t";
	for(auto iter = v2.begin(); iter!= v2.end(); iter++){
		std::cout << *iter << ",";
	}
	std::cout << "\t";
	for(auto iter = retTemp.first.begin(); iter!= retTemp.first.end(); iter++){
		std::cout << *iter << ",";
	}
	std::cout << "\t";
	for(auto iter = retTemp.second.begin(); iter!= retTemp.second.end(); iter++){
		std::cout << *iter << ",";
	}
	std::cout << std::endl;
	std::cout << std::to_string(retTemp.first.size() == 2) << std::endl;

	*/

	if(retTemp.first.size() == 2)
		return retTemp;
	else
		return ret;
}

// Find the symmetric difference of two vectors
std::vector<unsigned> utils::symmetricDiff(std::vector<unsigned> v1, std::vector<unsigned> v2, bool isSorted){
	std::vector<unsigned> ret;
	std::vector<unsigned> retTemp;

	if(v1 == v2)
		return ret;

	if(!isSorted){
		sort(v1.begin(), v1.end());
		sort(v2.begin(), v2.end());
	}
	set_symmetric_difference(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(retTemp));

	/*for(auto iter = v1.begin(); iter!= v1.end(); iter++){
		std::cout << *iter << ",";
	}
	std::cout << "\t";
	for(auto iter = v2.begin(); iter!= v2.end(); iter++){
		std::cout << *iter << ",";
	}
	std::cout << "\t";
	for(auto iter = retTemp.begin(); iter!= retTemp.end(); iter++){
		std::cout << *iter << ",";
	}
	std::cout << std::endl;*/

	return retTemp;
}

// Find the symmetric difference of two vectors
std::set<unsigned> utils::symmetricDiff(std::set<unsigned> v1, std::set<unsigned> v2, bool isSorted){
	std::set<unsigned> ret;
	std::set<unsigned> retTemp;

	if(v1 == v2)
		return ret;

	//if(!isSorted){
	//	sort(v1.begin(), v1.end());
	//	sort(v2.begin(), v2.end());
//	}
	set_symmetric_difference(v1.begin(), v1.end(), v2.begin(), v2.end(),  std::inserter(retTemp, retTemp.begin()));

	/*for(auto iter = v1.begin(); iter!= v1.end(); iter++){
		std::cout << *iter << ",";
	}
	std::cout << "\t";
	for(auto iter = v2.begin(); iter!= v2.end(); iter++){
		std::cout << *iter << ",";
	}
	std::cout << "\t";
	for(auto iter = retTemp.begin(); iter!= retTemp.end(); iter++){
		std::cout << *iter << ",";
	}
	std::cout << std::endl;*/

	return retTemp;
}

//Iteratively build subsets (faces) of the simplex set
std::vector<std::set<unsigned>> utils::getSubsets(std::set<unsigned> set, int dim){
	std::vector<std::set<unsigned>> subset;
	std::set<unsigned> empty;
	subset.push_back(empty);

	//For each set in the
	for(auto i = set.begin(); i!= set.end(); i++){
		std::vector<std::set<unsigned>> subsetTemp = subset;
		unsigned entry = *i;

		for (unsigned j = 0; j < subsetTemp.size(); j++){
			subsetTemp[j].insert(entry);
		}

		unsigned z = 0;
		for (auto j = subsetTemp.begin(); j != subsetTemp.end(); j++){
			subset.push_back(*j);

		}
	}

	std::vector<std::set<unsigned>> retSubset;

	for(std::set<unsigned> z : subset){
		if(z.size() == dim)
			retSubset.push_back(z);
	}
	return retSubset;
}

// Find the union of two vectors
std::set<unsigned> utils::setUnion(std::set<unsigned> v1, std::set<unsigned> v2){
	std::set<unsigned> retTemp;

	if(v1 == v2) return v1;

	set_union(v1.begin(), v1.end(), v2.begin(), v2.end(), std::inserter(retTemp, retTemp.begin()));

	return retTemp;
}



// Find the union of two vectors
std::vector<unsigned> utils::setUnion(std::vector<unsigned> v1, std::vector<unsigned> v2, bool isSorted){
	std::vector<unsigned> ret;
	std::vector<unsigned> retTemp;

	if(v1 == v2)
		return ret;

	if(!isSorted){
		sort(v1.begin(), v1.end());
		sort(v2.begin(), v2.end());
	}

	set_union(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(retTemp));

	/*for(auto iter = v1.begin(); iter!= v1.end(); iter++){
		std::cout << *iter << ",";
	}
	std::cout << "\t";
	for(auto iter = v2.begin(); iter!= v2.end(); iter++){
		std::cout << *iter << ",";
	}
	std::cout << "\t";
	for(auto iter = retTemp.begin(); iter!= retTemp.end(); iter++){
		std::cout << *iter << ",";
	}
	std::cout << std::endl;*/
	return retTemp;
}

void utils::writeLog(std::string module, std::string message){
	if(debug == "1" || debug == "true"){
		std::cout << "[" << module << "]:\t" << message << std::endl;
	} else {
		writeFile("[" + module + "]:\t" + message);
	}
	return;
}

void utils::writeDebug(std::string module, std::string message){
	if(debug == "0" || debug == "false"){
		return;
	} else {
		std::cout << "[DEBUG]\t[" << module << "]:\t" << message << std::endl;
		writeFile("[DEBUG]\t[" + module + "]:\t" + message);
	}

	return;
}

void utils::writeFile(std::string fullMessage){
	std::ofstream outfile;
	outfile.open(outputFile+"_debug.txt", std::ios_base::app);
	outfile << fullMessage << "\n";

	return;
}


bool utils::sortBySecond(const std::pair<std::set<unsigned>, double> &a, const std::pair<std::set<unsigned>, double> &b){
	if(a.second == b.second){ //If the simplices have the same weight, sort by reverse lexicographic order for fastPersistence
		auto itA = a.first.rbegin(), itB = b.first.rbegin();
		while(itA != a.first.rend()){
			if(*itA != *itB) return *itA>*itB;
			++itA; ++itB;
		}
		return false;
	} else{
		return a.second<b.second;
	}
}

//Iteratively build subsets (faces) of the simplex set
std::vector<std::set<unsigned>> utils::getSubsets(std::set<unsigned> set){
	std::vector<std::set<unsigned>> subset;
	std::set<unsigned> empty;
	subset.push_back(empty);

	//For each set in the
	for(auto i = set.begin(); i!= set.end(); i++){
		std::vector<std::set<unsigned>> subsetTemp = subset;
		unsigned entry = *i;

		for (unsigned j = 0; j < subsetTemp.size(); j++){
			subsetTemp[j].insert(entry);
		}

		unsigned z = 0;
		for (auto j = subsetTemp.begin(); j != subsetTemp.end(); j++){
			subset.push_back(*j);

		}
	}

	std::vector<std::set<unsigned>> retSubset;

	for(std::set<unsigned> z : subset){
		if(z.size() == set.size() - 1)
			retSubset.push_back(z);
	}
	return retSubset;
}

//Iteratively build subsets (faces) of the simplex set
std::vector<std::vector<unsigned>> utils::getSubsets(std::vector<unsigned> set){
	std::vector<std::vector<unsigned>> subset;
	std::vector<unsigned> empty;
	subset.push_back(empty);

	//For each set in the
	for(auto i = set.begin(); i!= set.end(); i++){
		std::vector<std::vector<unsigned>> subsetTemp = subset;
		unsigned entry = *i;

		for (unsigned j = 0; j < subsetTemp.size(); j++){
			subsetTemp[j].push_back(entry);
		}

		unsigned z = 0;
		for (auto j = subsetTemp.begin(); j != subsetTemp.end(); j++){
			subset.push_back(*j);

		}
	}

	std::vector<std::vector<unsigned>> retSubset;

	for(std::vector<unsigned> z : subset){
		if(z.size() == set.size() - 1)
			retSubset.push_back(z);
	}
	return retSubset;
}


std::vector<double> utils::nearestNeighbors(std::vector<double>& point, std::vector<std::vector<double>>& pointcloud){
	//based on random projection, x is current point being examined, n is number of centroids/facilities
	utils ut;
	std::vector<double> retVal;

	//Get sq distances for each point
	for(auto currentPoint : pointcloud){
		retVal.push_back(ut.vectors_distance(point, currentPoint));
	}

	return retVal;

}


std::vector<std::vector<double>> utils::deserialize(std::vector<double> serialData, unsigned dim){

	//First check if the vector size matches the dimension
	if(serialData.size() % dim != 0){
		std::cout << "Error occurred when deserializing data: invalid size" << std::endl;
		return {};
	}

	//Deduce the number of vectors
	auto n = serialData.size() / dim;
	auto begin = serialData.begin();

	std::vector<std::vector<double>> ret(n, std::vector<double>(dim));

	for(unsigned i = 0; i < n; i++){
		for(unsigned j = 0; j < dim; j++){
			ret[i][j] = (serialData)[(i*dim) + j];
		}
	}
	return ret;
}


std::vector<double> utils::serialize(std::vector<std::vector<double>> &origData){

	//Make sure we have data to serialize
	if(origData.size() == 0){
		std::cout << "Error occurred when serializing data: empty vector argument" << std::endl;
		return {};
	}

	auto n = origData.size();
	auto d = origData[0].size();

	//Preallocate our vector to prevent resizing
	std::vector<double> ret(n*d);

	//Copy element by element
	for(unsigned i = 0; i < n; i++){
		for(unsigned k = 0; k < d; k++){
			ret[(i*d)+k] = origData[i][k];
		}
	}

	return ret;
}


std::pair<std::vector<std::vector<double>>,std::vector<double>> utils::calculateBetaCentersandRadius(std::vector<unsigned> dsimplex ,std::vector<std::vector<double>> &inputData,std::vector<std::vector<double>>* distMatrix, double beta)
{
	std::vector<std::vector<double>> betacenters;
	std::vector<double> betaradii;
	bool intersectionCircle= false;
	bool intersectionLune = false;
	if(beta <0)
		exit(0);
	else if(beta==0)
		return std::make_pair(betacenters,betaradii);
    	else if(beta <1)
		intersectionCircle = true;
	else
		intersectionLune = true;

		if(beta<1)
			beta=1/beta;
    		std::set<unsigned> simplex(dsimplex.begin(),dsimplex.end());
	     	std::vector<double> circumCenter;
	   	if(simplex.size()>2)
			circumCenter = utils::circumCenter(simplex,inputData);
		else if(simplex.size()==2){
 			auto first = simplex.begin();
			std::vector<double> R;
			std::vector<double> A = inputData[*first];
			std::advance(first, 1);
			std::vector<double> B = inputData[*first];
			std::transform(A.begin(), A.end(), B.begin(), std::back_inserter(R),[](double e1,double e2){return ((e1+e2)/2);});
			circumCenter = R;
   		}
                
		double circumRadius;
		if(simplex.size()>2)
			circumRadius = utils::circumRadius(simplex,distMatrix);
		else
			circumRadius = pow((*distMatrix)[dsimplex[0]][dsimplex[1]]/2,2);
		bool first = true;
		
		std::vector<size_t> neighbors;
		std::vector<std::vector<size_t>> neighborsCircleIntersection;
		for (auto x : simplex){
			double expr1,expr2,expr3;
			std::vector<unsigned> face1;
			face1 = dsimplex;
			face1.erase(std::remove(face1.begin(),face1.end(),x),face1.end());
			std::set<unsigned> face(face1.begin(),face1.end());
			std::vector<double> faceCC ;
			if(face.size()>2)
				faceCC = utils::circumCenter(face,inputData);
			else if(face.size()==2){
				auto first = face.begin();
				std::vector<double> fR;
			}	
			double faceRadius;
   			if(face.size()>2)
   				faceRadius = utils::circumRadius(face,distMatrix);
   			else
				faceRadius = pow((*distMatrix)[face1[0]][face1[1]]/2,2);
			auto result = utils::nullSpaceOfMatrix(face,inputData,faceCC,sqrt(faceRadius));            
			std::vector<double> hpcoff = result.first;
			std::vector<double> betaCenter;
			double betaRadius;
             		std::vector<std::vector<double>> betaCenters ;
			bool sameside = false;
			if(intersectionCircle && beta >2){
				if(beta < 3){
					double ratio = sqrt(circumRadius)/sqrt(faceRadius);
		                        betaRadius = sqrt(faceRadius) + (beta-2)*(sqrt(circumRadius)-sqrt(faceRadius));
					betaCenters = utils::betaCentersCalculation(hpcoff, 1+(beta-2)*(ratio-1), sqrt(faceRadius),faceCC);
				}
				else{
					betaCenters = utils::betaCentersCalculation(hpcoff, beta-1, sqrt(faceRadius),faceCC);
		                 	betaRadius = sqrt(faceRadius)*(beta-1);
                                }
                  		expr1=0;
		  		expr2=0;
                                expr3 =0;
			      	for(unsigned i =0;i<hpcoff.size();i++){
				      	expr1 += hpcoff[i]*circumCenter[i];
					expr2 += hpcoff[i]*betaCenters[0][i];
					expr3 += hpcoff[i]*inputData[x][i];
		      		}
		  		expr1--;
		  		expr2--;
				expr3--;		
				if((expr1> 0 && expr3>0)&&expr2>0||(expr1<0&&expr3<0)&&expr2<0){
					sameside=true;
					betaCenter = betaCenters[1];
					
				}
				else if((expr1>0&&expr3>0)||(expr1<0&&expr3<0)){
					sameside=true;
					if(expr2>0&&expr1>0)
						betaCenter = betaCenters[1];
					else
						betaCenter = betaCenters[0];
				}else{
					sameside=false;
					if(expr2>0&&expr1>0)
						betaCenter = betaCenters[0];
					else
						betaCenter = betaCenters[1];
				}
                        }
			else{
                  		expr1=0;
		  		expr3=0;
			      	for(unsigned i =0;i<hpcoff.size();i++){
				      	expr1 += hpcoff[i]*circumCenter[i];
					expr3 += hpcoff[i]*inputData[x][i];
		      		}
		  		expr1--;
				expr3--;		
				if((expr1>0&&expr3>0)||(expr1<0&&expr3<0))
					sameside=true;
				else
					sameside=false;
					for(unsigned y =0 ;y< inputData[0].size();y++)
					if(sameside){
						if(intersectionCircle)
							betaCenter.push_back((2-beta)*circumCenter[y] + (beta-1)*faceCC[y]);
						else
							betaCenter.push_back(beta*circumCenter[y] - (beta-1)*faceCC[y]);
					}else{
						if(intersectionCircle)
						{
							betaCenter.push_back(beta*circumCenter[y] - (beta-1)*faceCC[y]);
						}
						else
							betaCenter.push_back((2-beta)*circumCenter[y] + (beta-1)*faceCC[y]);
           				}
			}
	   		if(!intersectionCircle || beta<=2)
		     	        betaRadius = utils::vectors_distance(betaCenter,inputData[face1[0]]);		
			
			betacenters.push_back(betaCenter);
			betaradii.push_back(betaRadius);
		}
		return std::make_pair(betacenters,betaradii);
}
