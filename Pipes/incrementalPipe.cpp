#include "incrementalPipe.hpp"
#include "alphaComplex.hpp"
#include "utils.hpp"
#include <fstream>
#include <omp.h>

double determinantOfMatrix(std::vector<std::vector<double>> mat, int n){
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

	return (det);}
std::vector<std::vector<double>> matrixMultiplication(std::vector<std::vector<double>> matA, std::vector<std::vector<double>> matB){
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
	return mat;}
std::vector<std::vector<double>> inverseOfMatrix(std::vector<std::vector<double>> mat, int n){
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
  return matinv;}
double circumRadius(std::set<unsigned> simplex,std::vector<std::vector<double>>* distMatrix){
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

	return -(determinantOfMatrix(matA,simplex.size())/(2*determinantOfMatrix(matACap,simplex.size()+1)));}
std::vector<double> circumCenter(std::set<unsigned> simplex,std::vector<std::vector<double>> inputData){
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



		return circumCenter;}

unsigned expand_d_minus_1_simplex()
{
return 0;
}

// basePipe constructor
template <typename nodeType>
incrementalPipe<nodeType>::incrementalPipe()
{
  this->pipeType = "incrementalPipe";
  return;
}

// runPipe -> Run the configured functions of this pipeline segment
template <typename nodeType>
void incrementalPipe<nodeType>::runPipe(pipePacket<nodeType> &inData)
{
	this->inputData=inData.inputData;
	unsigned dim = 17;
	this->search_space=std::vector<unsigned>(inData.inputData.size());
	std::vector<std::vector<unsigned>> dsimplexes;
	std::vector<std::vector<unsigned>> inner_dsimplexes_shell={{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16}};
	std::vector<std::vector<unsigned>> outer_dsimplexes_shell;
	std::vector<unsigned> simplex;
	while(search_space.size()!=0)
	{
		outer_dsimplexes_shell.clear();
		for(int i=0;i<1;i++)
		{
			for(unsigned j=0;j<dim;j++)
			{
				simplex=inner_dsimplexes_shell[i];
				simplex.erase(simplex.begin()+j);
				for(auto j: simplex){
					std::cout<<j<<" ";
				}
				std::cout<<std::endl;
			}
		}
		break;
	}
	return;
}

// configPipe -> configure the function settings of this pipeline segment
template <typename nodeType>
bool incrementalPipe<nodeType>::configPipe(std::map<std::string, std::string> &configMap)
{
  std::string strDebug;

  auto pipe = configMap.find("debug");
  if (pipe != configMap.end())
  {
    this->debug = std::atoi(configMap["debug"].c_str());
    strDebug = configMap["debug"];
  }
  pipe = configMap.find("outputFile");
  if (pipe != configMap.end())
    this->outputFile = configMap["outputFile"].c_str();

  this->ut = utils(strDebug, this->outputFile);

  this->configured = true;
  this->ut.writeDebug("delaunayPipe", "Configured with parameters { eps: " + configMap["epsilon"] + " , debug: " + strDebug + ", outputFile: " + this->outputFile + " }");

  return true;
}
// outputData -> used for tracking each stage of the pipeline's data output without runtime
template <typename nodeType>
void incrementalPipe<nodeType>::outputData(pipePacket<nodeType> &inData)
{
  std::ofstream file;
  file.open("output/" + this->pipeType + "_output.csv");
  // code to print the data

  file.close();
  return;
}

template class incrementalPipe<simplexNode>;
template class incrementalPipe<alphaNode>;
template class incrementalPipe<witnessNode>;