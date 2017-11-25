//Streams for input/output
#include <fstream>
#include <sstream>
#include <iostream>

//Datatypes
#include <string>
#include <utility> //std::pair
#include <tuple> //std::tie

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
				if (s!="edge"){
					std::cerr << "The input file does not obey dimacs format. Problem type is not 'edge'.";
					infile.close();
					return std::pair<int, ED::Graph> (EXIT_FAILURE, ED::Graph{0});
				}
				else {
					iss >> n >> m;
					goto initialization; //get out of the while loop
				}
				}//end case
			default : std::cerr << "The input file does not obey dimacs format.";
				infile.close();
				return std::pair<int, ED::Graph> (EXIT_FAILURE, ED::Graph{0});
		}
	}

	initialization: ;
	//Initialize the graph
	if (n==0){
		std::cerr << "There is a problem with the input file, e.g. no 'p' line or empty graph.";
		infile.close();
		return std::pair<int, ED::Graph> (EXIT_FAILURE, ED::Graph{0});
	}
	ED::Graph graph{(ED::NodeId) n};

	//Add the edges to the graph
	int edge_counter = 0;
	while(std::getline(infile, line)){
		iss.str (line);
		iss.clear();
		iss >> c;
		switch (c) {
			case 'c': break;
			case 'e': {
				int v, w;
				iss >> v >> w;
				graph.add_edge(ED::from_dimacs_id(v), ED::from_dimacs_id(w));
				edge_counter++;
				break;
				}
			default : std::cerr << "The input file does not obey dimacs format.";
				infile.close();
				return std::pair<int, ED::Graph> (EXIT_FAILURE, ED::Graph{0});
		}
		
	}

	infile.close();
	if (edge_counter != m){
		std::cerr << "The number of edges stated in the dimacs file is incorrect.";
		return std::pair<int, ED::Graph> (EXIT_FAILURE, graph);
	}
	return std::pair<int, ED::Graph> (EXIT_SUCCESS, graph);
}

}//end of namespace ED

//ToDo: To be removed
/*
int main(int argc, char* argv[])
{
	if (argc == 2){
		int success;
		ED::Graph graph(0);
		std::tie(success, graph) = ED::read_graph(argv[1]);
//		std::cout << "success: " + std::to_string(success) + "\n";
//		std::cout << graph;
		return 0;
	}
	else{
		std::cout << "Wrong number of arguments. One argument filename expected.";
		return 1;
	}
}
*/
