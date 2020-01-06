#include <iostream>
#include <vector>
#include "simplexBase.hpp"



// TEST simplexBase Functions
void t_simp_functions(std::string& log){
	simplexBase* testComplex = new simplexBase();
	std::string failLog = "";
	std::vector<double> testValue = {0.0, 1.0, 2.0};
	std::vector<std::vector<double>> testValueArray {{0.0, 1.0, 2.0},{2.0, 1.0, 0.0}};
	std::vector<unsigned> findValue = {0};
	std::set<unsigned> findValueSet = {0,1};
	
	//Insert values into uninitialized complex
	//	RET: void
	try{ testComplex->insert(testValue); }
	catch(const std::exception){ failLog += "simplexBase insert failed\n"; }
	
	//Get uninitialized size
	//	RET: double (-1)
	if(testComplex->getSize() != -1) { failLog += "simplexBase getSize failed\n"; }
	
	//Insert iterative into uninitialized complex
	//	RET: void
	if(testComplex->insertIterative(testValue, testValueArray)) { failLog += "simplexBase insertIterative failed\n"; }
	
	//Delete iterative from uninitialized complex
	//	RET: void
	try{ testComplex->deleteIterative(0); }
	catch(const std::exception){ failLog += "simplexBase deleteIterative failed \n"; }
	
	//Find vector unsigned in uninitialized complex
	//	RET: bool
	if(testComplex->find(findValue)) { failLog += "simplexBase find vector failed\n"; }
	
	//Find set unsigned in uninitialized complex
	//	RET: bool
	if(testComplex->find(findValueSet)) { failLog += "simplexBase find set failed\n"; }
	
	//Get uninitialized simplex count
	//	RET: int (-1)
	if(testComplex->simplexCount() != -1) { failLog += "simplexBase simplexCount failed\n"; }
	
	//Get uninitialized vertex count
	//	RET: int (-1)
	if(testComplex->vertexCount() != -1) { failLog += "simplexBase vertexCount failed\n"; }
	
	//Get dimensional edges from uninitialized complex
	//	RET: std::vector<std::vector<unsigned>>
	if(testComplex->getDimEdges(2, 5.0).size() != 0) { failLog += "simplexBase getDimEdges failed\n"; }
	
	//Get all edges from uninitialized complex
	//	RET: std::vector<std::vector<std::pair<std::set<unsigned>, double>>>
	if(testComplex->getAllEdges(5.0).size() != 0) { failLog += "simplexBase getAllEdges failed\n"; }
	
	//Get indexed edges from uninitialized complex
	//	RET: std::vector<std::vector<graphEntry>>
	if(testComplex->getIndexEdges(1).size() != 0) { failLog += "simplexBase getIndexEdges failed\n"; }
	
	//Expand the dimensions of the uninitialized complex
	//	RET: void
	try{ testComplex->expandDimensions(2); }
	catch(const std::exception){ failLog += "simplexBase expandDimensions failed\n"; }
	
	//Attempt to reduce the uninitialized complex
	//	RET: void
	try{ testComplex->reduceComplex(); }
	catch(const std::exception){ failLog += "simplexBase reduceComplex failed\n"; }
	
	//Attempt to trigger the stream evaluator for the uninitialized complex
	//	RET: bool
	testComplex->streamEvaluator(testValue, testValueArray);
	
	
	//Output log status to calling function
	if(failLog.size() > 0){
		log += "FAILED: simplexBase Test Functions---------------------------\n" + failLog;	
	} else {
		 log += "PASSED: simplexBase Test Functions---------------------------\n";
	}

	return;
}


// TEST simplexArrayList Functions
void t_simpArrayList_functions(){
	std::map<std::string,std::string> config;
	simplexBase* testComplex = new simplexBase();
	testComplex = testComplex->newSimplex("simplexArrayList", config);
	
}



int main (int, char**){
	std::string failLog;
	t_simp_functions(failLog);


	std::cout << std::endl << failLog << std::endl;
}
