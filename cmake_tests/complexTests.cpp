#include <iostream>
#include <vector>
#include "simplexBase.hpp"



// TEST simplexBase Functions
void t_simp_functions(std::string &log){
	/*simplexBase* testComplex = new simplexBase();
	std::string failLog = "";
	std::vector<double> testValue = {0.0, 1.0, 2.0};
	std::vector<std::vector<double>> testValueArray {{0.0, 1.0, 2.0},{2.0, 1.0, 0.0}};
	std::vector<unsigned> findValue = {0};
	std::set<unsigned> findValueSet = {0,1};
	// int keyToBeDeleted = 0;
	// int indexToBeDeleted = 0;
	// std::vector<double> distsFromCurrVec = {0.0, 2.83};

	std::cout << "Beginning simplexBase tests" << std::endl;
	//Insert values into uninitialized complex
	//	RET: void
	std::cout << "\tTesting insertion into uninitialized complex" << std::endl;
	try{ testComplex->insert(); }
	catch(const std::exception){ failLog += "simplexBase insert failed\n"; }

	//Get uninitialized size
	//	RET: double (-1)
	std::cout << "\tTesting get uninitialized size" << std::endl;
	if(testComplex->getSize() != -1) { failLog += "simplexBase getSize failed\n"; }

	//Insert iterative into uninitialized complex
	//	RET: void
	std::cout << "\tTesting insert iterative into uninitialized complex" << std::endl;
	if(testComplex->insertIterative(testValue, testValueArray)) { failLog += "simplexBase insertIterative failed\n"; }

	//Delete iterative from uninitialized complex
	//	RET: void
	std::cout << "\tTesting delete iterative from uninitialized complex" << std::endl;
	try{ testComplex->deleteIterative(0); }
	catch(const std::exception){ failLog += "simplexBase deleteIterative failed \n"; }

	//Find vector unsigned in uninitialized complex
	//	RET: bool
	std::cout << "\tTesting find vector unsigned in uninitialized complex" << std::endl;
	if(testComplex->find(findValue)) { failLog += "simplexBase find vector failed\n"; }

	//Find set unsigned in uninitialized complex
	//	RET: bool
	std::cout << "\tTesting find set unsigned in uninitialized complex" << std::endl;
	if(testComplex->find(findValueSet)) { failLog += "simplexBase find set failed\n"; }

	//Get uninitialized simplex count
	//	RET: int (-1)
	std::cout << "\tTesting get uninitialized simplex count" << std::endl;
	if(testComplex->simplexCount() != -1) { failLog += "simplexBase simplexCount failed\n"; }

	//Get uninitialized vertex count
	//	RET: int (-1)
	std::cout << "\tTesting get uninitialized vertex count" << std::endl;
	if(testComplex->vertexCount() != -1) { failLog += "simplexBase vertexCount failed\n"; }

	//Get dimensional edges from uninitialized complex
	//	RET: std::vector<std::vector<unsigned>>
	std::cout << "\tTesting get dimensional edges from uninitialized complex" << std::endl;
	if(testComplex->getDimEdges(2).size() != 0) { failLog += "simplexBase getDimEdges failed\n"; }

	//Get all edges from uninitialized complex
	//	RET: std::vector<std::vector<std::pair<std::set<unsigned>, double>>>
	std::cout << "\tTesting get all edges from uninitialized complex" << std::endl;
	if(testComplex->getAllEdges().size() != 0) { failLog += "simplexBase getAllEdges failed\n"; }

	//Expand the dimensions of the uninitialized complex
	//	RET: void
	std::cout << "\tTesting expand the dimensions of the uninitialized complex" << std::endl;
	//try{ testComplex->expandDimensions(2); }
	//catch(const std::exception){ failLog += "simplexBase expandDimensions failed\n"; }

	//Attempt to reduce the uninitialized complex
	//	RET: void
	std::cout << "\tTesting attempt to reduce the uninitialized complex" << std::endl;
	try{ testComplex->reduceComplex(); }
	catch(const std::exception){ failLog += "simplexBase reduceComplex failed\n"; }

	//Attempt to trigger the stream evaluator for the uninitialized complex
	//	RET: bool
	std::cout << "\tTesting attempt to trigger the stream evaluator for the uninitialized complex" << std::endl;
	testComplex->streamEvaluator(testValue, testValueArray);

	//Output log status to calling function
	if(failLog.size() > 0){
		log += "FAILED: simplexBase Test Functions---------------------------\n" + failLog;
	} else {
		 log += "PASSED: simplexBase Test Functions---------------------------\n";
	}
	*/
	std::cout << "\tFinished simplexBase tests" << std::endl;
	return;
}


// Test simplexTree insertion, removal, etc. functions for sliding window
void t_simp_tree_functions(std::string &log){
	
	/*
	std::map<std::string, std::string> config;
	std::string failLog = "";
	std::vector<std::vector<double>> testValueArray = {{0.0, 1.0, 2.0},{2.0, 1.0, 0.0}, {1.0, 1.0, 2.0}, \
													 {1.1, 1.1, 1.2},{0.0, 0.4, 1.0}, {1.5, 1.5, 0.0}	};

	std::vector<std::vector<double>> distMatrix = {};
	simplexBase* testComplex = simplexBase::newSimplex("simplexTree", config);

	
	std::cout << "Beginning simplexTree tests" << std::endl;
	
	//Insert values from testValueArray into the complex to build the distance matrix
	for(unsigned i = 0; i < testValueArray.size(); i++){
		//Insert iterative into simplexTree
		//	RET: void
		std::cout << "\tTesting insertIterative for simplexTree, index " << i << std::endl;
		if(testComplex->insertIterative(testValueArray[i], distMatrix)) { failLog += "simplexTree insertIterative failed\n"; }
	}

	//Attempt to delete values from the complex
	std::cout << "\tTesting deleteIterative for simplexTree, 0" << std::endl;
	testComplex->deleteIterative(0);
	std::cout << "\tTesting deleteIterative for simplexTree, 5" << std::endl;
	testComplex->deleteIterative(5);
	std::cout << "\tTesting deleteIterative for simplexTree, 3" << std::endl;
	testComplex->deleteIterative(3);

	//Output log status to calling function
	if(failLog.size() > 0){
		log += "FAILED: simplexTree Test Functions---------------------------\n" + failLog;
	} else {
		 log += "PASSED: simplexTree Test Functions---------------------------\n";
	}

	std::cout << "\tFinished simplexTree tests" << std::endl;
}



// TEST simplexArrayList Functions
void t_simp_base_functions(std::string &log, std::string type){
	std::map<std::string,std::string> config;
	std::string failLog = "";
	std::vector<double> testValue = {0.0, 1.0, 2.0};
	std::vector<std::vector<double>> testValueArray = {{0.0, 1.0, 2.0},{2.0, 1.0, 0.0}, {1.0, 1.0, 2.0}, \
													 {1.1, 1.1, 1.2},{0.0, 0.4, 1.0}, {1.5, 1.5, 0.0}	};
	std::vector<std::vector<double>> emptyValueArray = {};

	std::vector<unsigned> findValue = {0};
	std::set<unsigned> findValueSet = {0};

	simplexBase* testComplex = simplexBase::newSimplex(type, config);
	testComplex->setDistanceMatrix(&testValueArray);
	
	std::cout << "Beginning simplexBase functions for " << type << std::endl;

	//Insert values to initialize the complex
	//	RET: void
	try{
		std::cout << "\t" << type << "\tTesting insertion of values to initialize complex" << std::endl;
		for(auto vector : testValueArray)
			testComplex->insert();
	}
	catch(const std::exception){ failLog += type + " insert failed\n"; }

	//Get complex size
	//	RET: double (testValueArray.size())
	std::cout << "\t" << type << "\tTesting get complex size" << std::endl;
	if(testComplex->getSize() <= 0) { failLog += type + " getSize failed : " + std::to_string(testComplex->getSize()) + "\n"; }

	//Insert iterative into complex
	//	RET: void
	std::cout << "\t" << type << "\tTesting insert iterative into complex" << std::endl;
	if(testComplex->insertIterative(testValue, emptyValueArray)) { failLog += type + " insertIterative failed\n"; }

	//Delete iterative from complex
	//	RET: void
	std::cout << "\t" << type << "\tTesting delete iterative from complex" << std::endl;
	try{ testComplex->deleteIterative(0); }
	catch(const std::exception){ failLog += type + " deleteIterative failed \n"; }

	//Find vector unsigned in complex
	//	RET: bool
	std::cout << "\t" << type << "\tTesting find vector unsigned in complex" << std::endl;
	if(testComplex->find(findValue)) { failLog += type + " find vector failed\n"; }

	//Find set unsigned in complex
	//	RET: bool
	std::cout << "\t" << type << "\tTesting find set unsigned in complex" << std::endl;
	if(testComplex->find(findValueSet)) { failLog += type + " find set failed\n"; }

	//Get simplex count
	//	RET: int (-1)
	std::cout << "\t" << type << "\tTesting get simplex count" << std::endl;
	if(testComplex->simplexCount() != -1) { failLog += type + " simplexCount failed\n"; }

	//Get vertex count
	//	RET: int (-1)
	std::cout << "\t" << type << "\tTesting get vertex count" << std::endl;
	if(testComplex->vertexCount() != -1) { failLog += type + " vertexCount failed\n"; }

	//Get dimensional edges from complex
	//	RET: std::vector<std::vector<unsigned>>
	std::cout << "\t" << type << "\tTesting get dimensional edges from complex" << std::endl;
	std::cout << "\t" << testComplex->simplexList.size() << std::endl;
	if(testComplex->getDimEdges(2).size() != 0) { failLog += type + " getDimEdges failed\n"; }

	//Get all edges from complex
	//	RET: std::vector<std::vector<std::pair<std::set<unsigned>, double>>>
	std::cout << "\t" << type << "\tTesting get all edges from complex" << std::endl;
	if(testComplex->getAllEdges().size() != 0) { failLog += type + " getAllEdges failed\n"; }

	//Expand the dimensions of the complex
	//	RET: void
	std::cout << "\t" << type << "\tTesting expand the dimensions of the complex" << std::endl;
	//try{ testComplex->expandDimensions(2); }
	//catch(const std::exception){ failLog += type + " expandDimensions failed\n"; }

	//Attempt to reduce the complex
	//	RET: void
	std::cout << "\t" << type << "\tTesting attempt to reduce the complex" << std::endl;
	//try{ testComplex->reduceComplex(); }
	//catch(const std::exception){ failLog += type + " reduceComplex failed\n"; }

	//Attempt to trigger the stream evaluator for the complex
	//	RET: bool
	std::cout << "\t" << type << "\tTesting attempt to trigger the stream evaluator for the complex" << std::endl;
	//testComplex->streamEvaluator(testValue, testValueArray);

	//Output log status to calling function
	if(failLog.size() > 0){
		log += "FAILED: " + type + " Base Test Functions---------------------------\n" + failLog;
	} else {
		 log += "PASSED: " + type + " Base Test Functions---------------------------\n";
	}
	
	std::cout << "\tFinished simplexBase functions for " << type << std::endl;
	*/
}

// TEST simplexArrayList Functions
void t_simp_empty_functions(std::string &log, std::string type){
	/*
	std::map<std::string,std::string> config;
	std::string failLog = "";
	std::vector<double> testValue = {0.0, 1.0, 2.0};
	std::vector<std::vector<double>> testValueArray {{0.0, 1.0, 2.0},{2.0, 1.0, 0.0}};
	std::vector<unsigned> findValue = {0};
	std::set<unsigned> findValueSet = {0};
	
	std::cout << "Beginning simplexBase empty functions for " << type << std::endl;

	simplexBase* testComplex = simplexBase::newSimplex(type, config);
	//Don't insert anything into the complex and test our functions

	//Get empty complex size
	//	RET: double (0)
	std::cout << "\t" << type << "\tTesting get size from empty complex" << std::endl;
	if(testComplex->getSize() != 0) { failLog += type + " getSize failed : " + std::to_string(testComplex->getSize()) + "\n"; }

	//Delete iterative from empty complex
	//	RET: void
	//std::cout << "\t" << type << "\tTesting delete iterative from empty complex" << std::endl;
	//try{ testComplex->deleteIterative(0,0); }
	//catch(const std::exception){ failLog += type + " deleteIterative failed \n"; }

	//Find vector unsigned in empty complex
	//	RET: bool (false)
	std::cout << "\t" << type << "\tTesting find vector unsigned from empty complex" << std::endl;
	if(testComplex->find(findValue)) { failLog += type + " find vector failed\n"; }

	//Find set unsigned in empty complex
	//	RET: bool (false)
	std::cout << "\t" << type << "\tTesting find set unsigned from empty complex" << std::endl;
	if(testComplex->find(findValueSet)) { failLog += type + " find set failed\n"; }

	//Get empty simplex count
	//	RET: int (0)
	std::cout << "\t" << type << "\tTesting get simplex count from empty complex" << std::endl;
	if(testComplex->simplexCount() != 0) { failLog += type + " simplexCount failed\n"; }

	//Get empty vertex count
	//	RET: int (0)
	std::cout << "\t" << type << "\tTesting get vertex count from empty complex" << std::endl;
	if(testComplex->vertexCount() != 0) { failLog += type + " vertexCount failed\n"; }

	//Get dimensional edges from empty complex
	//	RET: std::vector<std::vector<unsigned>>
	std::cout << "\t" << type << "\tTesting get dim edges from empty complex" << std::endl;
	if(testComplex->getDimEdges(2).size() != 0) { failLog += type + " getDimEdges failed\n"; }

	//Get all edges from empty complex
	//	RET: std::vector<std::vector<std::pair<std::set<unsigned>, double>>>
	std::cout << "\t" << type << "\tTesting get all edges from empty complex" << std::endl;
	if(testComplex->getAllEdges().size() != 0) { failLog += type + " getAllEdges failed\n"; }

	//Expand the dimensions of the empty complex
	//	RET: void
	std::cout << "\t" << type << "\tTesting expand dimensions from empty complex" << std::endl;
	//try{ testComplex->expandDimensions(2); }
	//catch(const std::exception){ failLog += type + " expandDimensions failed\n"; }

	//Attempt to reduce the empty complex
	//	RET: void
	std::cout << "\t" << type << "\tTesting reduce complex from empty complex" << std::endl;
	try{ testComplex->reduceComplex(); }
	catch(const std::exception){ failLog += type + " reduceComplex failed\n"; }

	//Attempt to trigger the stream evaluator for the empty complex
	//	RET: bool
	std::cout << "\t" << type << "\tTesting stream evaluator from empty complex" << std::endl;
	testComplex->streamEvaluator(testValue, testValueArray);

	//Output log status to calling function
	if(failLog.size() > 0){
		log += "FAILED: " + type + " Empty Test Functions---------------------------\n" + failLog;
	} else {
		 log += "PASSED: " + type + " Empty Test Functions---------------------------\n";
	}
	
	std::cout << "\tFinished simplexBase empty functions for " << type << std::endl;
	*/
}

int main (int, char**){
	std::string log;
	t_simp_functions(log);

	for(std::string type : {"simplexArrayList","simplexTree"}){
		try{t_simp_empty_functions(log, type);}
		catch(const std::exception){log += "FAILED: " + type + " Empty Test Functions---------------------------\n";}

		//try{t_simp_base_functions(log, type);}
		//catch(const std::exception){log += "FAILED: " + type + " Base Test Functions---------------------------\n";}

	}

	std::cout << std::endl << log << std::endl;
}

