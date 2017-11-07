#include <iostream>
#include <cstdlib>

#include "graph.cpp"
#include "graph.hpp"

//return graph consisting of path with \c num_nodes many vertices
/*static ED::Graph create_path(ED::NodeId num_nodes)
{
   ED::Graph result{num_nodes};

   for (ED::NodeId node_id = 0; node_id + 1 < num_nodes; ++node_id)
   {
      result.add_edge(node_id, node_id + 1);
   }
   //will be moved
   return result;
}*/

void agument_M(ED::Graph M, std::vector <ED::NodeId> path) {
	int i=0;
	for (unsigned int j=0; j != path.size()-2; ++j) {
		if(i==0) {
			M.add_edge(path[j],path[j+1]);
			i=1;
		}
		if(i==1) {
			M.delete_edge(path[j], path[j+1]);
			i=0;
		}
	}
}

bool DFSI(ED::Graph const & T, ED::NodeId v, ED::NodeId w, std::vector <ED::NodeId> & path, std::vector <bool> & visited){
	visited[v]=1;
	auto const & node = T.node(v);
	path.push_back(v);
	for (auto const & neighbor_id : node.neighbors()) {
		if(neighbor_id==w) return 1;
		else if(visited[neighbor_id]==0 /*& visited[R[neighbor_id]]==0*/){
			if(DFSI(T, neighbor_id, w, path, visited)) return 1;
			path.pop_back();
		}
		/*else if(visited[neighbor_id]==0 & R[neighbor_Id]==R[v]) {
			if(DFSI(T, neighbor_id, w, path, visited)) return 1;
			path.pop_back();
		}*/
	} 
	return 0;
}

//Depth first search
std::vector <ED::NodeId> DFS(ED::Graph const & T, const int size,  ED::NodeId v, ED::NodeId w) {
	std::vector <ED::NodeId> path;
	std::vector <bool> visited(size, 0);
	DFSI(T,v,w,path,visited);
	return path;
}


//checks if v has a neighbour w with w\notin Odd(T), if not it returns v
ED::NodeId check_neighbours_of_node(ED::NodeId v, ED::Graph const & G, int T[]) {
	auto const & node = G.node(v);
	for (auto const & neighbor_id : node.neighbors()) {
		if(T[neighbor_id]!=2) {
			return neighbor_id;
		}
	}
	return v;
}

//checks if v is part of T
bool check_node_in_T(ED::NodeId v, int T[]) {
	if(T[v]==0){
		return 0;
	}
	else return 1;
}

//checks if v is a vertrex of T or covered by the matching M
bool check_node_M_exposed(ED::NodeId v, ED::Graph M) {
	auto const & node = M.node(v);
	if (&node.neighbors().size!=0) {
		return 0;
	}
	return 1;
}

//checks if M is a perfect matching
bool check_M_perfect(ED::Graph M){
	int i=0;
	for (int node_id = 0; node_id < (int)M.num_nodes(); ++node_id) {
		auto const & node = M.node(node_id);
		if (&node.neighbors().size!=0) {
			i=1;
		}
		if(i==0) return 0;
		i=0;
	}
	return 1;
}

ED::Graph perfect_Matching(ED::Graph G) {
	ED::Graph M{G.num_nodes()};
	
	if(G.num_nodes()==0 || G.num_nodes()==1){ 
		std::cout<< "G is empty or has only one vertex";
		return M;
	}
	//we choose two structures for the tree
	//first is a graph Tree with num_nodes() verticies and, at the beginning, no edges
	//second is a vector T with num_nodes() entries
	//we set T[x]=0 if x\in Tree, T[x]=1 if x\in Even(T), T[x]=2 if x\in Odd(T)
	
	ED::Graph Tree{G.num_nodes()};
	std::vector <int> T(G.num_nodes());
	
	ED::NodeId start=0;  // root of T
	T[start]=1;
	
	//to simulate the shrinking we use partition classes
	std::vector <ED::NodeId> Partition_class(G.num_nodes());
	for(unsigned int i=0; i<G.num_nodes(); i++){
		Partition_class[i]=i;		//at the beginning each vertex has its own partition class 
	} 
	
	return G;
}


int main()
{	
	/*ED::Graph g1{4};
	int T[g1.num_nodes()];
   // create and output a path with argc many vertices
    //ED::Graph const graph = create_path(static_cast<ED::NodeId>(argc));
   
	g1.add_edge(1,2);
	g1.add_edge(1,3);
	auto const & node = g1.node(0);
	for (auto const & neighbor_id : node.neighbors()) {
		std::cout<<neighbor_id;
	}
	int i;
	for(i=0;i<5;++i) {
		std::cout<<i;
	}

   std::cout << graph;
   std::cout << g1;
   return EXIT_SUCCESS;
   */
   return 1;
}
