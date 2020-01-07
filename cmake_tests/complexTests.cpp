#include <iostream>
#include <vector>
#include "simplexBase.hpp"



// TEST simplexBase Functions
void t_simp_functions(std::string &log){
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
void t_simp_base_functions(std::string &log, std::string type){
	std::map<std::string,std::string> config;
	simplexBase* testComplex = new simplexBase();
	std::string failLog = "";
	std::vector<double> testValue = {0.0, 1.0, 2.0};
	std::vector<std::vector<double>> testValueArray {{0.0, 1.0, 2.0},{2.0, 1.0, 0.0}, {1.0, 1.0, 2.0}, \
													 {1.1, 1.1, 1.2},{0.0, 0.4, 1.0}, {1.5, 1.5, 0.0}	};
	std::vector<unsigned> findValue = {0};
	std::set<unsigned> findValueSet = {0};
	
	testComplex = testComplex->newSimplex(type, config);
	//Insert values to initialize the complex
	//	RET: void
	try{
		for(auto vector : testValueArray)
			testComplex->insert(vector); 
	}
	catch(const std::exception){ failLog += type + " insert failed\n"; }
	
	//Get complex size
	//	RET: double (testValueArray.size())
	if(testComplex->getSize() <= 0) { failLog += type + " getSize failed : " + std::to_string(testComplex->getSize()) + "\n"; }
	
	//Insert iterative into uninitialized complex
	//	RET: void
	if(testComplex->insertIterative(testValue, testValueArray)) { failLog += type + " insertIterative failed\n"; }
	
	//Delete iterative from uninitialized complex
	//	RET: void
	try{ testComplex->deleteIterative(0); }
	catch(const std::exception){ failLog += type + " deleteIterative failed \n"; }
	
	//Find vector unsigned in uninitialized complex
	//	RET: bool
	if(testComplex->find(findValue)) { failLog += type + " find vector failed\n"; }
	
	//Find set unsigned in uninitialized complex
	//	RET: bool
	if(testComplex->find(findValueSet)) { failLog += type + " find set failed\n"; }
	
	//Get uninitialized simplex count
	//	RET: int (-1)
	if(testComplex->simplexCount() != -1) { failLog += type + " simplexCount failed\n"; }
	
	//Get uninitialized vertex count
	//	RET: int (-1)
	if(testComplex->vertexCount() != -1) { failLog += type + " vertexCount failed\n"; }
	
	//Get dimensional edges from uninitialized complex
	//	RET: std::vector<std::vector<unsigned>>
	if(testComplex->getDimEdges(2, 5.0).size() != 0) { failLog += type + " getDimEdges failed\n"; }
	
	//Get all edges from uninitialized complex
	//	RET: std::vector<std::vector<std::pair<std::set<unsigned>, double>>>
	if(testComplex->getAllEdges(5.0).size() != 0) { failLog += type + " getAllEdges failed\n"; }
	
	//Get indexed edges from uninitialized complex
	//	RET: std::vector<std::vector<graphEntry>>
	if(testComplex->getIndexEdges(1).size() != 0) { failLog += type + " getIndexEdges failed\n"; }
	
	//Expand the dimensions of the uninitialized complex
	//	RET: void
	try{ testComplex->expandDimensions(2); }
	catch(const std::exception){ failLog += type + " expandDimensions failed\n"; }
	
	//Attempt to reduce the uninitialized complex
	//	RET: void
	try{ testComplex->reduceComplex(); }
	catch(const std::exception){ failLog += type + " reduceComplex failed\n"; }
	
	//Attempt to trigger the stream evaluator for the uninitialized complex
	//	RET: bool
	testComplex->streamEvaluator(testValue, testValueArray);
	
	
	//Output log status to calling function
	if(failLog.size() > 0){
		log += "FAILED: " + type + " Base Test Functions---------------------------\n" + failLog;	
	} else {
		 log += "PASSED: " + type + " Base Test Functions---------------------------\n";
	}
	
}


// TEST simplexArrayList Functions
void t_simp_empty_functions(std::string &log, std::string type){
	std::map<std::string,std::string> config;
	simplexBase* testComplex = new simplexBase();
	std::string failLog = "";
	std::vector<double> testValue = {0.0, 1.0, 2.0};
	std::vector<std::vector<double>> testValueArray {{0.0, 1.0, 2.0},{2.0, 1.0, 0.0}};
	std::vector<unsigned> findValue = {0};
	std::set<unsigned> findValueSet = {0};
	
	testComplex = testComplex->newSimplex(type, config);
	//Don't insert anything into the complex and test our functions
	
	//Get empty complex size
	//	RET: double (0)
	if(testComplex->getSize() != 0) { failLog += type + " getSize failed : " + std::to_string(testComplex->getSize()) + "\n"; }
	
	//Delete iterative from empty complex
	//	RET: void
	try{ testComplex->deleteIterative(0); }
	catch(const std::exception){ failLog += type + " deleteIterative failed \n"; }
	
	//Find vector unsigned in empty complex
	//	RET: bool (false)
	if(testComplex->find(findValue)) { failLog += type + " find vector failed\n"; }
	//Find set unsigned in empty complex
	//	RET: bool (false)
	if(testComplex->find(findValueSet)) { failLog += type + " find set failed\n"; }
	//Get empty simplex count
	//	RET: int (0)
	if(testComplex->simplexCount() != 0) { failLog += type + " simplexCount failed\n"; }
	
	//Get empty vertex count
	//	RET: int (0)
	if(testComplex->vertexCount() != 0) { failLog += type + " vertexCount failed\n"; }
	
	//Get dimensional edges from empty complex
	//	RET: std::vector<std::vector<unsigned>>
	if(testComplex->getDimEdges(2, 5.0).size() != 0) { failLog += type + " getDimEdges failed\n"; }
	
	//Get all edges from empty complex
	//	RET: std::vector<std::vector<std::pair<std::set<unsigned>, double>>>
	if(testComplex->getAllEdges(5.0).size() != 0) { failLog += type + " getAllEdges failed\n"; }
	
	//Get indexed edges from empty complex
	//	RET: std::vector<std::vector<graphEntry>>
	if(testComplex->getIndexEdges(1).size() != 0) { failLog += type + " getIndexEdges failed\n"; }
	
	//Expand the dimensions of the empty complex
	//	RET: void
	try{ testComplex->expandDimensions(2); }
	catch(const std::exception){ failLog += type + " expandDimensions failed\n"; }
	
	//Attempt to reduce the empty complex
	//	RET: void
	try{ testComplex->reduceComplex(); }
	catch(const std::exception){ failLog += type + " reduceComplex failed\n"; }
	
	//Attempt to trigger the stream evaluator for the empty complex
	//	RET: bool
	testComplex->streamEvaluator(testValue, testValueArray);
	
	//Output log status to calling function
	if(failLog.size() > 0){
		log += "FAILED: " + type + " Empty Test Functions---------------------------\n" + failLog;	
	} else {
		 log += "PASSED: " + type + " Empty Test Functions---------------------------\n";
	}
}

int main (int, char**){
	std::string log;
	t_simp_functions(log);
	
	for(std::string type : {"simplexArrayList"}){ //,"simplexTree","indSimplexTree"}){ TODO
		try{t_simp_empty_functions(log, type);}
		catch(const std::exception){log += "FAILED: " + type + " Empty Test Functions---------------------------\n";}	
		
		try{t_simp_base_functions(log, type);}
		catch(const std::exception){log += "FAILED: " + type + " Base Test Functions---------------------------\n";}	
		
	}
	
	std::cout << std::endl << log << std::endl;
}
