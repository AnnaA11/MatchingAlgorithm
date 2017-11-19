//Streams for input/output
#include <fstream>
#include <sstream>
#include <iostream>

//Datatypes
#include <string>
#include <utility> //std::pair

//Graph class
#include "graph.cpp"
#include "graph.hpp"

namespace ED
{

/*
Reads in a graph from the file filename, which has to be in DIMACS format
@return: pair of int success and ED::Graph

*/
std::pair<int, ED::Graph> read_graph(std::string filename)
{
	//Open the file filename
	std::ifstream infile(filename);
	std::string line;
	
	if (!infile){
		std::cerr << "Could not read file "+filename;
		return std::pair <int, ED::Graph> (EXIT_FAILURE, ED::Graph {0});
	}

	std::istringstream iss;
	char c;
	int n = 0;
	int m = 0;

	//Find the "problem" line and extract n and m
	while(std::getline(infile, line)){
		iss.str (line);
		iss.clear();
		iss >> c;
		switch (c) {
			case 'c': break;
			case 'p': {
				std::string s;
				iss >> s;
				if (s=="edge"){
					iss >> n >> m;
					break;
				}
				else {
					std::cerr << "The input file does not obey dimacs format.";
					infile.close();
					return std::pair<int, ED::Graph> (EXIT_FAILURE, ED::Graph{0});
				}
				}//end case
			default : std::cerr << "The input file does not obey dimacs format.";
				infile.close();
				return std::pair<int, ED::Graph> (EXIT_FAILURE, ED::Graph{0});
		}
	}

	//Initialize the graph
	if (n==0){
		std::cerr << "There is a problem with the input file, e.g. no 'p' line or empty graph.";
		infile.close();
		return std::pair<int, ED::Graph> (EXIT_FAILURE, ED::Graph{0});
	}
	ED::Graph graph{(ED::NodeId) n};

	//Add the edges to the graph
	while(std::getline(infile, line)){
		iss.str (line);
		iss.clear();
		iss >> c;
		switch (c) {
			case 'c': break;
			case 'e': {
				char v,w;
				iss >> v >> w;
				graph.add_edge(ED::from_dimacs_id(v), ED::from_dimacs_id(w));
				break;
				}
			default : std::cerr << "The input file does not obey dimacs format.";
				infile.close();
				return std::pair<int, ED::Graph> (EXIT_FAILURE, ED::Graph{0});
		}
		
	}

	infile.close();
	return std::pair<int, ED::Graph> (EXIT_SUCCESS, graph);
}

}//end of namespace ED

//ToDo: To be removed
int main()
{ return 0; }
