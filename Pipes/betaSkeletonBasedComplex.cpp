/*
 * betaSkeletonBasedComplexPipe hpp + cpp extend the basePipe class for calculating the
 * beta Skeleton Based Complex generation for data input
 *
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
#include "betaSkeletonBasedComplex.hpp"
#include "alphaComplex.hpp"
#include "qhullPipe.hpp"
#include "utils.hpp"
#include "readInput.hpp"

// basePipe constructor
template <typename nodeType>
betaSkeletonBasedComplex<nodeType>::betaSkeletonBasedComplex()
{
	this->pipeType = "betaSkeletonBasedComplex";
	return;
}

// runPipe -> Run the configured functions of this pipeline segment
template <typename nodeType>
void betaSkeletonBasedComplex<nodeType>::runPipe(pipePacket<nodeType> &inData)
{
	// Generate Beta Skeleton Based Complex
	// Temporarily commenting this out - need to check inData.complex type
	//		If type is the graph-based simplexArrayList (inherited) then
	//			cast to gAL and run non-virtual function:

	//	((graphArrayList*)inData.complex)->graphInducedComplex(dim,inData.inputData,beta);
	std::vector<std::vector<unsigned>> dsimplexmesh;
	kdTree tree(inData.inputData, inData.inputData.size()); // KDTree for efficient nearest neighbor search
	int dim = inData.inputData[0].size();
	double distanceSum = 0;
	int count = 0;
	// Computing the average point cloud distance
	for (int x = 0; x < inData.inputData.size(); x++)
		for (int y = x + 1; y < inData.inputData.size(); y++)
		{
			distanceSum += (*((alphaComplex<nodeType> *)inData.complex)->distMatrix)[x][y];
			count = count + 1;
		}
	double averageDistance = distanceSum / count;
	count = 0;
	std::vector<std::pair<int, int>> neighborsepsilon;

	// Enumeration within epsilon ball for beta <1.
	if (this->betaMesh == "null.csv")
	{
		if (this->beta < 1)
		{
			for (unsigned index = 0; index < inData.inputData.size(); index++)
			{
				std::vector<size_t> neighbors = tree.neighborhoodIndices(inData.inputData[index], this->epsilon); // All neighbors in epsilon-ball
				int n = neighbors.size();
				neighborsepsilon.push_back(std::make_pair(n, index));
			}
			sort(neighborsepsilon.begin(), neighborsepsilon.end());
			std::vector<int> toremove;
			for (unsigned indexold = 0; indexold < inData.inputData.size(); indexold++)
			{
				toremove.push_back(neighborsepsilon[indexold].second);
				unsigned index = neighborsepsilon[indexold].second;
				std::vector<size_t> neighbors = tree.neighborhoodIndices(inData.inputData[index], this->epsilon); // All neighbors in epsilon-ball
				int n = neighbors.size();
				std::sort(toremove.begin(), toremove.end());
				std::sort(neighbors.begin(), neighbors.end());
				std::vector<int> difference;
				std::set_difference(neighbors.begin(), neighbors.end(), toremove.begin(), toremove.end(), std::back_inserter(difference));
				n = difference.size();
				if (n >= dim)
				{
					std::vector<unsigned> dsimplex(dim);
					std::vector<unsigned> dsimplexIndexed;
					std::vector<unsigned>::iterator first = dsimplex.begin(), last = dsimplex.end();
					std::generate(first, last, UniqueNumber);
					dsimplexIndexed.push_back(index);
					for (int i = 0; i < dim; i++)
						dsimplexIndexed.push_back(difference[dsimplex[i] - 1]);
					// for each enumerated simplex evaluate for beta selecton criterea
					if (checkInsertDsimplex(dsimplexIndexed, inData, this->beta, averageDistance, tree))
					{
						std::sort(dsimplexIndexed.begin(), dsimplexIndexed.end());
						dsimplexmesh.push_back(dsimplexIndexed);
						count++;
					}
					while ((*first) != n - dim + 1)
					{
						std::vector<unsigned>::iterator mt = last;
						while (*(--mt) == n - (last - mt) + 1)
							;
						(*mt)++;
						while (++mt != last)
							*mt = *(mt - 1) + 1;
						std::vector<unsigned> dsimplexIndexed1;
						dsimplexIndexed1.push_back(index);
						for (int i = 0; i < dim; i++)
						{
							dsimplexIndexed1.push_back(difference[dsimplex[i] - 1]);
							//   std::cout<<difference[dsimplex[i]-1]<<" ";
						}
						// for each enumerated simplex evaluate for beta selecton criterea
						if (checkInsertDsimplex(dsimplexIndexed1, inData, this->beta, averageDistance, tree))
						{
							std::sort(dsimplexIndexed1.begin(), dsimplexIndexed1.end());
							dsimplexmesh.push_back(dsimplexIndexed1);
							count++;
						}
					}
				}
			}
		}
		// Compute Delaunay Triangulation for beta >=1.
		else
		{
			Qhull qh;
			kdTree tree(inData.inputData, inData.inputData.size()); // KDTree for efficient nearest neighbor search
			std::vector<double> sdata;
			// serializing all the data
			for (auto a : inData.inputData)
				for (auto b : a)
					sdata.push_back(b);
			PointCoordinates *pts = new PointCoordinates(qh, inData.inputData[0].size(), "UCI Data Sets");
			pts->append(sdata);
			qh.runQhull(pts->comment().c_str(), pts->dimension(), pts->count(), &*pts->coordinates(), "d o");
			std::vector<std::vector<int>> dsimplexes1 = qdelaunay_o(qh);
			std::vector<std::vector<int>> surfacefacets = qconvex_o(qh);
			for (auto x : surfacefacets)
			{
				for (auto y : x)
					std::cout << y << " ";
				std::cout << "\n";
			}
			std::vector<std::vector<unsigned>> dsimplexes;
			for (auto x : dsimplexes1)
			{
				std::vector<unsigned> temp;
				for (auto y : x)
					temp.push_back(y);
				dsimplexes.push_back(temp);
			}
			for (auto x : dsimplexes)
				// for each enumerated simplex evaluate for beta selecton criterea
				if (checkInsertDsimplex(x, inData, this->beta, averageDistance, tree))
				{
					std::sort(x.begin(), x.end());
					dsimplexmesh.push_back(x);
				}
		}
		std::ofstream out("latestbetamesh.csv");
		for (auto x : dsimplexmesh)
		{
			int i = 0;
			for (auto y : x)
			{
				i++;
				out << y;
				if (i != x.size())
					out << ",";
			}
			out << '\n';
		}
	}
	else if (betaMesh != "null.csv")
	{
		auto rs = readInput();
		std::vector<std::vector<double>> betaMesh1 = rs.readCSV(this->betaMesh);
		std::vector<std::vector<unsigned>> betaMesh;
		for (auto x : betaMesh1)
		{
			std::vector<unsigned> temp;
			for (auto y : x)
				temp.push_back(y);
			betaMesh.push_back(temp);
		}
		for (auto x : betaMesh)
			// for each enumerated simplex evaluate for beta selecton criterea
			if (checkInsertDsimplex(x, inData, this->beta, averageDistance, tree))
			{
				std::sort(x.begin(), x.end());
				dsimplexmesh.push_back(x);
			}
		std::ofstream out("latestbetamesh.csv");
		for (auto x : dsimplexmesh)
		{
			int i = 0;
			for (auto y : x)
			{
				i++;
				out << y;
				if (i != x.size())
					out << ",";
			}
			out << '\n';
		}
	}

	//	((alphaComplex<nodeType>*)inData.complex)->buildBetaComplex(dsimplexmesh,inData.inputData.size(),inData.inputData,this->beta,this->betaMode);
	//	((alphaComplex<nodeType>*)inData.complex)->buildAlphaComplex(dsimplexmesh,inData.inputData.size(),inData.inputData);
	//	((alphaComplex<nodeType>*)inData.complex)->buildFilteration(dsimplexmesh,inData.inputData.size(),inData.inputData,this->beta);
	//        ((alphaComplex<nodeType>*)inData.complex)->buildBetaComplexFilteration(dsimplexmesh, inData.inputData.size(),inData.inputData, tree);
	std::vector<std::vector<bool>> incidenceMatrix(inData.inputData.size(), std::vector<bool>(inData.inputData.size(), 0));
	int countd1 = 0;
	for (auto x : dsimplexmesh)
	{
		for (int i = 0; i < x.size(); i++)
			for (int j = i + 1; j < x.size(); j++)
			{
				int origin = x[i];
				int destination = x[j];
				if (incidenceMatrix[origin][destination] != true)
					countd1++;
				incidenceMatrix[origin][destination] = true;
			}
	}

	// std::cout<<"Edges ="<<countd1<<" ";
	inData.incidenceMatrix = incidenceMatrix;
	std::ofstream file("PHdSphereDimensionWiseMeshSize.txt", std::ios_base::app);
	file << this->betaMode << "," << inData.inputData.size() << "," << inData.inputData[0].size() << "," << this->beta << "," << dsimplexmesh.size() << std::endl;
	file.close();
	this->ut.writeDebug("betaSkeletonBasedComplex Pipe", "\tbetaSkeletonBasedComplex Size: ");
	return;
}
template <typename nodeType>

std::vector<std::vector<int>> betaSkeletonBasedComplex<nodeType>::qconvex_o(const Qhull &qhull)
{
	int dim = qhull.hullDimension();
	int numfacets = qhull.facetList().count();
	int totneighbors = numfacets * dim; /* incorrect for non-simplicial facets, see qh_countfacets */
	std::cout << "dim " << dim << "\n"
			  << "qhull size" << qhull.points().size() << " number of facets " << numfacets << "totalneighbors/2 " << totneighbors / 2 << "\n";
	std::vector<std::vector<double>> points;
	// for(QhullPoint point : qhull.points())
	for (QhullPoints::ConstIterator i = qhull.points().begin(); i != qhull.points().end(); ++i)
	{
		QhullPoint point = *i;
		points.push_back(point.toStdVector());
	}
	for (auto x : points)
	{
		for (auto y : x)
			std::cout << y << " ";
		std::cout << "\n";
	}
	QhullFacetList facets = qhull.facetList();
	std::vector<std::vector<int>> facetVertices;
	// for(QhullFacet f : facets)
	QhullFacetListIterator j(facets);

	while (j.hasNext())
	{
		QhullFacet f = j.next();
		std::vector<int> vertices;
		if (!f.isGood())
		{
			// ignore facet
		}
		else if (!f.isTopOrient() && f.isSimplicial())
		{ /* orient the vertices like option 'o' */
			QhullVertexSet vs = f.vertices();
			vertices.push_back(vs[1].point().id());
			vertices.push_back(vs[0].point().id());
			for (int i = 2; i < (int)vs.size(); ++i)
			{
				vertices.push_back(vs[i].point().id());
			}
			facetVertices.push_back(vertices);
		}
		else
		{ /* note: for non-simplicial facets, this code does not duplicate option 'o', see qh_facet3vertex and qh_printfacetNvertex_nonsimplicial */
			// for(QhullVertex vertex : f.vertices()){
			QhullVertexSetIterator k(f.vertices());
			while (k.hasNext())
			{
				QhullVertex vertex = k.next();
				QhullPoint p = vertex.point();
				vertices.push_back(p.id());
			}
			facetVertices.push_back(vertices);
		}
	}
	return facetVertices;
} // qconvex_o

template <typename nodeType>
std::vector<std::vector<int>> betaSkeletonBasedComplex<nodeType>::qdelaunay_o(const Qhull &qhull)
{
	int hullDimension = qhull.hullDimension();
	std::vector<std::vector<double>> inputSites;
	QhullPoints points = qhull.points();

	QhullPointsIterator j(points);
	while (j.hasNext())
	{
		QhullPoint point = j.next();
		inputSites.push_back(point.toStdVector());
	}
	QhullFacetList facets = qhull.facetList();
	int numFacets = facets.count();
	size_t numRidges = numFacets * hullDimension / 2;

	std::vector<std::vector<int>> regions;
	QhullFacetListIterator k(facets);
	while (k.hasNext())
	{
		QhullFacet f = k.next();
		std::vector<int> vertices;
		if (!f.isUpperDelaunay())
		{
			if (!f.isTopOrient() && f.isSimplicial())
			{
				QhullVertexSet vs = f.vertices();
				vertices.push_back(vs[1].point().id());
				vertices.push_back(vs[0].point().id());
				for (int i = 2; i < (int)vs.size(); ++i)
				{
					vertices.push_back(vs[i].point().id());
				}
			}
			else
			{
				QhullVertexSetIterator i(f.vertices());
				while (i.hasNext())
				{
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
template <>
unsigned betaSkeletonBasedComplex<alphaNode>::selectCenter(std::vector<double> hpcofffaces, std::vector<std::vector<double>> betaCenters, std::vector<double> otherPoint)
{
	double valuebetaCenter1 = 0, valuebetaCenter2 = 0, valueotherPoint = 0;
	for (unsigned i = 0; i < hpcofffaces.size(); i++)
	{
		valuebetaCenter1 += hpcofffaces[i] * betaCenters[0][i];
		valuebetaCenter2 += hpcofffaces[i] * betaCenters[1][i];
		valueotherPoint += hpcofffaces[i] * otherPoint[i];
	}
	if (valueotherPoint < 0)
	{
		if (valuebetaCenter1 < 0)
			return 0;
		else if (valuebetaCenter2 < 0)
			return 1;
	}
	else if (valueotherPoint > 0)
	{
		if (valuebetaCenter1 > 0)
			return 0;
		else if (valuebetaCenter2 > 0)
			return 1;
	}
	return 1;
}

template <>
bool betaSkeletonBasedComplex<alphaNode>::checkCC_Simplex_Inclusion(std::vector<unsigned> simplex, std::vector<std::vector<double>> inputData, std::vector<double> circumCenter)
{
	int i = 0;
	std::vector<std::vector<double>> matT;
	std::vector<std::vector<double>> matPv;
	for (auto x : simplex)
	{
		if (i != 1)
		{
			std::vector<double> tempmat;
			for (int j = 0; j < inputData[0].size(); j++)
				tempmat.push_back(inputData[x][1] - inputData[x][j]);
			matT.push_back(tempmat);
		}
		if (i == 1)
		{
			for (int j = 0; j < inputData[0].size(); j++)
			{
				std::vector<double> tempmat;
				tempmat.push_back(circumCenter[j] - inputData[x][j]);
				matPv.push_back(tempmat);
			}
		}
		i++;
	}
	std::vector<std::vector<double>> transposematT(matT.size(), std::vector<double>(matT[0].size(), 0));
	for (int i = 0; i < matT.size(); ++i)
		for (int j = 0; j < matT[0].size(); ++j)
			transposematT[j][i] = matT[i][j];

	std::vector<std::vector<double>> lambda = utils::matrixMultiplication(transposematT, matPv);
	bool outside = false;
	double sum = 0;
	for (auto x : lambda)
	{
		sum += x[0];
		if (x[0] < 0)
		{
			outside = true;
			break;
		}
	}
	if (sum > 1)
		outside = true;
	return outside;
}
template <>
int betaSkeletonBasedComplex<alphaNode>::getoppvertex(std::vector<unsigned> simplex, std::vector<std::vector<double>> inputData, std::vector<double> circumCenter)
{
	int oppvertex = -1;
	double mindist = 999999;
	for (auto x : simplex)
	{
		std::vector<unsigned> face1;
		face1 = simplex;
		face1.erase(std::remove(face1.begin(), face1.end(), x), face1.end());
		std::set<unsigned> face(face1.begin(), face1.end());
		std::vector<double> faceCC;
		if (face.size() > 2)
			faceCC = utils::circumCenter(face, inputData);
		else if (face.size() == 2)
		{
			auto first = face.begin();
			std::vector<double> fR;
			std::vector<double> fA = inputData[*first];
			std::advance(first, 1);
			std::vector<double> fB = inputData[*first];
			std::transform(fA.begin(), fA.end(), fB.begin(), std::back_inserter(fR), [](double e1, double e2)
						   { return ((e1 + e2) / 2); });
			faceCC = fR;
		}
		if (mindist > utils::vectors_distance(faceCC, circumCenter))
			oppvertex = x;
	}
	return oppvertex;
}
template <>
bool betaSkeletonBasedComplex<alphaNode>::checkInsertDsimplex(std::vector<unsigned> dsimplex, pipePacket<alphaNode> &inData, double beta, double averageDistance, kdTree tree)
{
	double maxEdge = 0;
	for (auto x : dsimplex)
		for (auto y : dsimplex)
			if (maxEdge < (*((alphaComplex<alphaNode> *)inData.complex)->distMatrix)[x][y])
				maxEdge = (*((alphaComplex<alphaNode> *)inData.complex)->distMatrix)[x][y];

	if (maxEdge > this->epsilon)
		return false;

	std::vector<size_t> neighborsfinalLune;
	std::vector<size_t> neighborsfinalCircle;

	bool intersectionCircle = false;
	bool intersectionLune = false;
	if (beta < 0)
		exit(0);
	else if (beta == 0)
		return true;
	else if (beta < 1)
		intersectionCircle = true;
	else if (beta > 1)
		intersectionLune = true;
	else if (beta == 1)
	{
		intersectionLune = true;
		intersectionCircle = true;
	}
	/*
		if(this->betaMode == "highDimCirclePrevious"){
			if(beta<1)
				beta=1/beta;
				std::set<unsigned> simplex(dsimplex.begin(),dsimplex.end());
				std::vector<double> circumCenter;
			if(simplex.size()>2)
				circumCenter = utils::circumCenter(simplex,inData.inputData);
			else if(simplex.size()==2){
				auto first = simplex.begin();
				std::vector<double> R;
				std::vector<double> A = inData.inputData[*first];
				std::advance(first, 1);
				std::vector<double> B = inData.inputData[*first];
				std::transform(A.begin(), A.end(), B.begin(), std::back_inserter(R),[](double e1,double e2){return ((e1+e2)/2);});
				circumCenter = R;
			}

			double circumRadius;
			if(simplex.size()>2)
				circumRadius = utils::circumRadius(simplex,((alphaComplex<alphaNode>*)inData.complex)->distMatrix);
			else
				circumRadius = pow((*((alphaComplex<alphaNode>*)inData.complex)->distMatrix)[dsimplex[0]][dsimplex[1]]/2,2);
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
					faceCC = utils::circumCenter(face,inData.inputData);
				else if(face.size()==2){
					auto first = face.begin();
					std::vector<double> fR;
					std::vector<double> fA = inData.inputData[*first];
							std::advance(first, 1);
					std::vector<double> fB = inData.inputData[*first];
					std::transform(fA.begin(), fA.end(), fB.begin(), std::back_inserter(fR),[](double e1,double e2){return ((e1+e2)/2);});
					faceCC = fR;
				}
				double faceRadius;
				if(face.size()>2)
					faceRadius = utils::circumRadius(face,((alphaComplex<alphaNode>*)inData.complex)->distMatrix);
				else
					faceRadius = pow((*((alphaComplex<alphaNode>*)inData.complex)->distMatrix)[face1[0]][face1[1]]/2,2);

				std::vector<double> hpcoff = utils::nullSpaceOfMatrix(face,inData.inputData,faceCC,sqrt(faceRadius));
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
								betaRadius = sqrt(circumRadius);
									}
							expr1=0;
					expr2=0;
									expr3 =0;
						for(unsigned i =0;i<hpcoff.size();i++){
							expr1 += hpcoff[i]*circumCenter[i];
						expr2 += hpcoff[i]*betaCenters[0][i];
						expr3 += hpcoff[i]*inData.inputData[x][i];
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
						expr3 += hpcoff[i]*inData.inputData[x][i];
						}
					expr1--;
					expr3--;
					if((expr1>0&&expr3>0)||(expr1<0&&expr3<0))
						sameside=true;
					else
						sameside=false;
						for(unsigned y =0 ;y< inData.inputData[0].size();y++)
						if(sameside){
							if(intersectionCircle)
								betaCenter.push_back((2-beta)*circumCenter[y] + (beta-1)*faceCC[y]);
							else
								betaCenter.push_back(beta*circumCenter[y] - (beta-1)*faceCC[y]);
						}else{
							if(intersectionCircle)
								betaCenter.push_back(beta*circumCenter[y] - (beta-1)*faceCC[y]);
							else
								betaCenter.push_back((2-beta)*circumCenter[y] + (beta-1)*faceCC[y]);
							}
				}
				if(!intersectionCircle || beta<=2)
							betaRadius = utils::vectors_distance(betaCenter,inData.inputData[face1[0]]);

				std::cout<<betaRadius<<" ";
				for(auto x:betaCenter)
					std::cout<<x<<" ";
				std::cout<<"\n";
		//		int pk;
		//		std::cin>>pk;
				std::vector<size_t> neighborsface = tree.neighborhoodIndices(betaCenter, betaRadius); //All neighbors in epsilon-ball


						neighborsface.erase(std::remove(neighborsface.begin(),neighborsface.end(),x),neighborsface.end());
					std::sort (neighborsface.begin(),neighborsface.end());
					std::sort (neighbors.begin(),neighbors.end());
				  //	neighborsCircleIntersection.push_back(neighborsface);
				for(auto x: neighborsface)
				std::cout<<x<<" ";
				std::cout<<"\n";
				if(!first){
						 if(intersectionCircle == true){
						std::vector<size_t> v(std::min(neighbors.size(),neighborsface.size()));
						std::vector<size_t>::iterator it;
						it=std::set_intersection (neighbors.begin(), neighbors.end(), neighborsface.begin(), neighborsface.end(), v.begin());
						v.resize(it-v.begin());
						neighbors = v;
					}
					 else{
						std::vector<size_t> v(neighbors.size() + neighborsface.size());
						std::vector<size_t>::iterator it;
						it=std::set_union (neighbors.begin(), neighbors.end(), neighborsface.begin(), neighborsface.end(), v.begin());
						v.resize(it-v.begin());
						neighbors = v;
					}
					}
					else{
						neighbors = neighborsface;
						first = false;
					}
			}
			std::vector<size_t> circumneighbors = tree.neighborhoodIndices(circumCenter, sqrt(circumRadius)); //All neighbors in epsilon-ball
			for(auto y :dsimplex)
				circumneighbors.erase(std::remove(circumneighbors.begin(),circumneighbors.end(),y),circumneighbors.end());
			std::sort (circumneighbors.begin(),circumneighbors.end());
			bool circleintersect = false;
			bool brloop = false;
			for(int i=0;i<dsimplex.size()&&!brloop;i++)
				for(int j=i+1;j<dsimplex.size()&&!brloop;j++){
					std::vector<size_t> v(std::min(neighborsCircleIntersection[i].size(),neighborsCircleIntersection[j].size()));
					std::vector<size_t>::iterator it;
					it=std::set_intersection (neighborsCircleIntersection[i].begin(), neighborsCircleIntersection[i].end(), neighborsCircleIntersection[j].begin(), neighborsCircleIntersection[j].end(), v.begin());
					v.resize(it-v.begin());
					if(intersectionCircle)
					{
					std::vector<size_t> newv(std::min(circumneighbors.size(),v.size()));
					std::vector<size_t>::iterator it2;
					it2=std::set_intersection (circumneighbors.begin(), circumneighbors.end(), v.begin(), v.end(), newv.begin());
					newv.resize(it2-newv.begin());
						if(newv.size()>0){
						circleintersect = true;
						brloop = true;
					}
					}
					else{
						if(v.size()>0){
							circleintersect = true;
							brloop=true;
						}
					}
				}

						 //   if(intersectionCircle){
		//			if(!circleintersect && neighbors.size()==0)
		//				return true;
		//			else
		//				return false;
		//		}
		//		else if(!circleintersect&&neighbors.size()==0)
				if(neighbors.size()==0)
					return true;
				else
					return false;
		}
		*/

	// Checking the beta celtic knot check for given beta value  Circle mode is generally for beta <1
	if (this->betaMode == "highDimCircle")
	{
		if (beta < 1)
			beta = 1 / beta;
		std::set<unsigned> simplex(dsimplex.begin(), dsimplex.end());
		std::vector<double> circumCenter;
		if (simplex.size() > 2)
			circumCenter = utils::circumCenter(simplex, inData.inputData);
		else if (simplex.size() == 2)
		{
			auto first = simplex.begin();
			std::vector<double> R;
			std::vector<double> A = inData.inputData[*first];
			std::advance(first, 1);
			std::vector<double> B = inData.inputData[*first];
			std::transform(A.begin(), A.end(), B.begin(), std::back_inserter(R), [](double e1, double e2)
						   { return ((e1 + e2) / 2); });
			circumCenter = R;
		}

		double circumRadius;
		if (simplex.size() > 2)
			circumRadius = utils::circumRadius(simplex, ((alphaComplex<alphaNode> *)inData.complex)->distMatrix);
		else
			circumRadius = pow((*((alphaComplex<alphaNode> *)inData.complex)->distMatrix)[dsimplex[0]][dsimplex[1]] / 2, 2);
		bool first = true;

		bool obtuse = false;
		obtuse = checkCC_Simplex_Inclusion(dsimplex, inData.inputData, circumCenter);
		int CCfacingfacet = -1;
		if (obtuse)
		{
			CCfacingfacet = getoppvertex(dsimplex, inData.inputData, circumCenter);
		}
		std::vector<size_t> neighbors;
		std::vector<std::vector<size_t>> neighborsCircleIntersection;
		for (auto x : simplex)
		{
			std::vector<unsigned> face1;
			face1 = dsimplex;
			face1.erase(std::remove(face1.begin(), face1.end(), x), face1.end());
			std::set<unsigned> face(face1.begin(), face1.end());
			std::vector<double> faceCC;
			if (face.size() > 2)
				faceCC = utils::circumCenter(face, inData.inputData);
			else if (face.size() == 2)
			{
				auto first = face.begin();
				std::vector<double> fR;
				std::vector<double> fA = inData.inputData[*first];
				std::advance(first, 1);
				std::vector<double> fB = inData.inputData[*first];
				std::transform(fA.begin(), fA.end(), fB.begin(), std::back_inserter(fR), [](double e1, double e2)
							   { return ((e1 + e2) / 2); });
				faceCC = fR;
			}
			double faceRadius;
			if (face.size() > 2)
				faceRadius = utils::circumRadius(face, ((alphaComplex<alphaNode> *)inData.complex)->distMatrix);
			else
				faceRadius = pow(((*((alphaComplex<alphaNode> *)inData.complex)->distMatrix)[face1[0]][face1[1]] / 2), 2);
			auto result = utils::nullSpaceOfMatrix(face, inData.inputData, faceCC, sqrt(faceRadius));
			std::vector<double> hpcoff = result.first;
			std::vector<std::vector<double>> refbetaCenters;
			bool sameside = false;
			double fixdistance = utils::vectors_distance(faceCC, circumCenter) + sqrt(circumRadius);
			//	std::cout<<fixdistance<<" "<<sqrt(faceRadius)<<"  "<<sqrt(circumRadius)<<"\n";
			double refbeta = sqrt((circumRadius / faceRadius) + 1);
			//	std::cout<<"\n  "<<refbeta<<"\nBefore squareroot "<<std::abs(((fixdistance*fixdistance)/(faceRadius))-1)<<"\n";
			std::vector<double> refpoint;
			refbetaCenters = utils::betaCentersCalculation(hpcoff, refbeta, sqrt(faceRadius), faceCC);

			//***************************************************Obtuse facet facing circumcenter*************************************************************************************
			// double expr1=0;
			// double expr3=0;
			//  	for(unsigned i =0;i<hpcoff.size();i++){
			//      	expr1 += hpcoff[i]*refbetaCenters[0][i];
			//	    expr3 += hpcoff[i]*inData.inputData[x][i];
			//}
			// expr1--;
			// expr3--;
			//**************************************************Obtuse facet facing circumcenter***************************************************************************************
			if (obtuse && x == CCfacingfacet)
				sameside = false;
			else
				sameside = true;
			if (sameside)
				refpoint = refbetaCenters[0];
			else
				refpoint = refbetaCenters[1];

			//		for(auto x:circumCenter)
			//			std::cout<<x<<" ";
			std::vector<double> betaCenter;
			//		std::cout<<"refpoint ";
			//		for(auto x:refpoint)
			//			std::cout<<x<<" ";
			for (unsigned y = 0; y < inData.inputData[0].size(); y++)
				betaCenter.push_back(beta * circumCenter[y] + (1 - beta) * refpoint[y]);
			double betaRadius = utils::vectors_distance(betaCenter, inData.inputData[face1[0]]);

			//		std::cout<<"\n"<<betaRadius<<" ";
			//		for(auto x:betaCenter)
			//			std::cout<<x<<" ";
			//		std::cout<<"\nVertex "<<x<<"\n Dsimplex";
			//         std::cout<<beta<<" ";
			//		for(auto x:dsimplex)
			//			std::cout<<x<<" ";
			// 	           int pk;
			//     	std::cin>>pk;
			std::vector<size_t> neighbors1faces1 = tree.neighborhoodIndices(betaCenter, betaRadius); // All neighbors in epsilon-ball
			for (auto t : dsimplex)
				neighbors1faces1.erase(std::remove(neighbors1faces1.begin(), neighbors1faces1.end(), t), neighbors1faces1.end());
			//	for(auto x: neighbors1faces1)
			//	std::cout<<x<<" ";
			//	std::cout<<"\n";
			if (!first)
			{
				std::sort(neighborsfinalCircle.begin(), neighborsfinalCircle.end());
				std::sort(neighbors1faces1.begin(), neighbors1faces1.end());
				if (intersectionCircle == true)
				{
					std::vector<size_t> v1(std::min(neighborsfinalCircle.size(), neighbors1faces1.size()));
					std::vector<size_t>::iterator it1;
					it1 = std::set_intersection(neighbors1faces1.begin(), neighbors1faces1.end(), neighborsfinalCircle.begin(), neighborsfinalCircle.end(), v1.begin());
					v1.resize(it1 - v1.begin());
					neighborsfinalCircle = v1;
				}
				else
				{
					std::vector<size_t> v1(std::max(neighborsfinalCircle.size(), neighbors1faces1.size()));
					std::vector<size_t>::iterator it1;
					it1 = std::set_union(neighbors1faces1.begin(), neighbors1faces1.end(), neighborsfinalCircle.begin(), neighborsfinalCircle.end(), v1.begin());
					v1.resize(it1 - v1.begin());
					neighborsfinalCircle = v1;
				}
			}
			else
				neighborsfinalCircle = neighbors1faces1;
			first = false;
		}
	}
	// Checking the beta celtic knot check for given beta value  lune mode is generally for beta >=1
	else if (this->betaMode == "highDimLune")
	{
		std::set<unsigned> simplex(dsimplex.begin(), dsimplex.end());
		std::vector<double> circumCenter;
		std::vector<double> circumCenterfaces;
		std::vector<double> circumCenterfaces1;
		if (simplex.size() > 2)
			circumCenter = utils::circumCenter(simplex, inData.inputData);
		else if (simplex.size() == 2)
		{
			auto first = simplex.begin();
			std::vector<double> R;
			std::vector<double> A = inData.inputData[*first];
			std::advance(first, 1);
			std::vector<double> B = inData.inputData[*first];
			std::transform(A.begin(), A.end(), B.begin(), std::back_inserter(R), [](double e1, double e2)
						   { return ((e1 + e2) / 2); });
			circumCenter = R;
		}
		double circumRadius;
		if (simplex.size() > 2)
			circumRadius = utils::circumRadius(simplex, ((alphaComplex<alphaNode> *)inData.complex)->distMatrix);
		else
			circumRadius = pow((*((alphaComplex<alphaNode> *)inData.complex)->distMatrix)[dsimplex[0]][dsimplex[1]] / 2, 2);
		bool first = true;
		for (auto x : simplex)
		{
			std::vector<double> betaCenter;
			for (unsigned y = 0; y < inData.inputData[0].size(); y++)
				betaCenter.push_back(beta * circumCenter[y] + (1 - beta) * inData.inputData[x][y]);
			double betaRadius = beta * sqrt(circumRadius);

			/*	std::cout<<betaRadius<<" ";
				for(auto x:betaCenter)
					std::cout<<x<<" ";
				std::cout<<"\n";
			//	std::cout<<betaRadius<<" ";
		//		for(auto x:betaCenter)
		//			std::cout<<x<<" ";
		  //              int pk;
		//		std::cin>>pk;
		* */
			std::vector<size_t> neighbors1faces1 = tree.neighborhoodIndices(betaCenter, betaRadius); // All neighbors in epsilon-ball
			neighbors1faces1.erase(std::remove(neighbors1faces1.begin(), neighbors1faces1.end(), x), neighbors1faces1.end());
			/*		for(auto x: neighbors1faces1)
					std::cout<<x<<" ";
					std::cout<<"\n";
				*/
			if (!first)
			{
				std::sort(neighborsfinalLune.begin(), neighborsfinalLune.end());
				std::sort(neighbors1faces1.begin(), neighbors1faces1.end());
				if (intersectionLune == true)
				{
					std::vector<size_t> v1(std::min(neighborsfinalLune.size(), neighbors1faces1.size()));
					std::vector<size_t>::iterator it1;
					it1 = std::set_intersection(neighbors1faces1.begin(), neighbors1faces1.end(), neighborsfinalLune.begin(), neighborsfinalLune.end(), v1.begin());
					v1.resize(it1 - v1.begin());
					neighborsfinalLune = v1;
				}
				else
				{
					std::vector<size_t> v1(std::max(neighborsfinalLune.size(), neighbors1faces1.size()));
					std::vector<size_t>::iterator it1;
					it1 = std::set_union(neighbors1faces1.begin(), neighbors1faces1.end(), neighborsfinalLune.begin(), neighborsfinalLune.end(), v1.begin());
					v1.resize(it1 - v1.begin());
					neighborsfinalLune = v1;
				}
			}
			else
				neighborsfinalLune = neighbors1faces1;
			first = false;
		}
	}
	// std::cout<<this->betaMode ;
	if (this->betaMode == "highDimLune" && neighborsfinalLune.size() == 0)
		return true;
	if (this->betaMode == "highDimCircle" && neighborsfinalCircle.size() == 0)
	{
		return true;
	}
	return false;
}

// configPipe -> configure the function settings of this pipeline segment
template <typename nodeType>
bool betaSkeletonBasedComplex<nodeType>::configPipe(std::map<std::string, std::string> &configMap)
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

	pipe = configMap.find("beta");
	if (pipe != configMap.end())
		this->beta = std::atof(configMap["beta"].c_str());

	pipe = configMap.find("betaMode");
	if (pipe != configMap.end())
		this->betaMode = configMap["betaMode"].c_str();

	pipe = configMap.find("epsilon");
	if (pipe != configMap.end())
		this->epsilon = std::atof(configMap["epsilon"].c_str());

	this->ut = utils(strDebug, this->outputFile);
	pipe = configMap.find("dimensions");
	if (pipe != configMap.end())
	{
		this->dim = std::atoi(configMap["dimensions"].c_str());
	}
	pipe = configMap.find("betaMesh");
	if (pipe != configMap.end())
		this->betaMesh = configMap["betaMesh"].c_str();

	pipe = configMap.find("epsilon");
	if (pipe != configMap.end())
		this->enclosingRadius = std::atof(configMap["epsilon"].c_str());
	else
		return false;

	this->configured = true;
	this->ut.writeDebug("betaSkeletonBasedComplex Pipe ", "Configured with parameters { eps: " + configMap["epsilon"] + configMap["beta"] + " , debug: " + strDebug + ", outputFile: " + this->outputFile + " }");

	return true;
}

// outputData -> used for tracking each stage of the pipeline's data output without runtime
template <typename nodeType>
void betaSkeletonBasedComplex<nodeType>::outputData(pipePacket<nodeType> &inData)
{
	// Output related to betaSkeletonBasedComplex
	return;
}
template <typename nodeType>
bool betaSkeletonBasedComplex<nodeType>::checkInsertDsimplex(std::vector<unsigned> dsimplex, pipePacket<nodeType> &inData, double beta, double averageDistance, kdTree tree)
{
	std::cout << "No function Defined";
	return false;
}

template class betaSkeletonBasedComplex<simplexNode>;
template class betaSkeletonBasedComplex<alphaNode>;
template class betaSkeletonBasedComplex<witnessNode>;
