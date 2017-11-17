#include <iostream>
#include <cstdlib>

#include "graph.cpp"
#include "graph.hpp"


void agument_M(ED::Graph & M, std::vector <ED::NodeId> const & path) {
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



bool path_M_augmenting(ED::Graph const & M, std::vector <ED::NodeId> const & path) {
	int i=0;
	for(int j=0; j< path.size(); j++){
		if(i==1){
			if(M.node(path.at(j)).neighbors().at(0)!=path.at(j+1)) return 0;
		}
		else {
			i=1-i;
		}
	}
	return 1;
}

// DFS only computes a path in T, it has not to be M augmenting in the partition classes
// so we compute in each partition class the M augmenting path
bool find_M_agumenting_path_inclass(ED::Graph const & T, ED::Graph const & M, std::vector <bool> & visited, std::vector <int> const & Partition, std::vector <ED::NodeId> & p, ED::		NodeId start, ED::NodeId end) {
	visited.at(start)=1;
	p.push_back(start);
	for(auto const & neighbor_id : (T.node(start)).neighbors()) {
		if(M.node(neighbor_id).neighbors().at(0)==end){
			p.push_back(neighbor_id);
			p.push_back(M.node(neighbor_id).neighbors().at(0));
			return 1;
		}
		else if(Partition[neighbor_id]==Partition[start] && visited.at(neighbor_id)!=0) {
			if(Partition[M.node(neighbor_id).neighbors().at(0)]==Partition[start]) {
				p.push_back(neighbor_id);
				if(find_M_agumenting_path_inclass(T,M,visited,Partition,p,M.node(neighbor_id).neighbors().at(0),end)==1) return 1;
				p.pop_back();
				p.pop_back();
			}
		}
	}
	return 0;
}

std::vector <ED::NodeId> find_M_agumenting_path(ED::Graph const & M, ED::Graph const & T, std::vector <int> const & Partition, std::vector <ED::NodeId> const & p) {
	int start=0;
	int end=0;
	int k=Partition[p.at(0)];
	std::vector <bool> visited(M.num_nodes(), 0);
	
	std::vector <ED::NodeId> path;
	path.push_back(p.at(0));
	
	for(int j=1; j<p.size(); j++) {
		if(k==Partition[p.at(j)]){
			end=j;
		}
		else {
			if(end!=start){
				find_M_agumenting_path_inclass(T, M, visited, Partition, path, p.at(start), p.at(end));
			}
			else path.push_back(p.at(start));
			k=Partition[p.at(j)];
			start=j;
			end=j;
		}
	}
	return path;
}

bool DFSI(ED::Graph const & T, ED::NodeId v, ED::NodeId w, std::vector <ED::NodeId> & path, std::vector <bool> & visited, bool i){
	visited[v]=1;
	auto const & node = T.node(v);
	path.push_back(v);
	for (auto const & neighbor_id : node.neighbors()) {
		if(neighbor_id==w) return 1;
		else if(visited[neighbor_id]==0 /*&& visited[R[neighbor_id]]==0*/){    // eventuell unnötig
			if(DFSI(T, neighbor_id, w, path, visited,1-i)) return 1;
			path.pop_back();
		}
		/*else if(visited[neighbor_id]==0 && R[neighbor_Id]==R[v]) {		// eventuell unnötig
			if(DFSI(T, neighbor_id, w, path, visited,1-i)) return 1;
			path.pop_back();
		}*/
	} 
	return 0;
}

//Depth first search
std::vector <ED::NodeId> DFS(ED::Graph const & T, const int size,  ED::NodeId v, ED::NodeId w) {
	std::vector <ED::NodeId> path;
	bool i=1;
	std::vector <bool> visited(size, 0);
	DFSI(T,v,w,path,visited, i);
	return path;
}


//checks if v has a neighbour w with w\notin Odd(T), if not it returns v
ED::NodeId check_neighbours_of_node(ED::NodeId v, ED::Graph const & G, std::vector <int> const & T, std::vector <int> const & Partition) {
	auto const & node = G.node(v);
	for (auto const & neighbor_id : node.neighbors()) {
		if(T[neighbor_id]!=2 && Partition[neighbor_id]!=Partition[v]) {
			return neighbor_id;
		}
	}
	return v;
}


//checks if v is covered by the matching M
bool check_node_M_exposed(ED::NodeId v, ED::Graph const & M) {
	auto const & node = M.node(v);
	if (node.neighbors().size()!=0) {
		return 0;
	}
	return 1;
}

bool find_M_exposed_node(ED::Graph const & M, ED::NodeId v){
	for(unsigned int i=0; i<M.num_nodes(); i++){
		if(&(M.node(i)).neighbors().size==0){
			v=i;
			return 1;
		}
	}
	return 0;
}

//checks if M is a perfect matching
bool check_M_perfect(ED::Graph const & M){
	for (int node_id = 0; node_id < (int)M.num_nodes(); ++node_id) {
		if (((M.node(node_id)).neighbors()).size()==0) {
			return 0;
		}
	}
	return 1;
}

//we choose three structures for the tree
//first is a graph Tree with num_nodes() verticies and, at the beginning, no edges
//second is an int vector T with num_nodes() entries
//we set T[x]=0 if x\in Tree, T[x]=1 if x\in Even(T), T[x]=2 if x\in Odd(T)ED::Graph Tree{G.num_nodes()};
//third is a vector which containe all even verticies in T

//grow tree corresponds to one iteration of the while loop in the perfect matching algorithm
int grow_tree (ED::Graph const & G, ED::Graph & M, ED::NodeId start, ED::Graph & Tree, std::vector <int> & T, std::vector <ED::NodeId> & T_even, 
				std::vector <int> & Partition, std::vector <std::vector <ED::NodeId>> & Partition_classes) {
	ED::NodeId x,y;
	
	for(int i=0; i<T_even.size(); i++){
		
		y=check_neighbours_of_node(T_even.at(i), G, T, Partition);
		
		if(y!=T_even.at(i)) {
			x=i;
			if(T[y]==0 && check_node_M_exposed(y,M)){
				agument_M(M,find_M_agumenting_path(M,Tree,Partition,DFS(Tree, Tree.num_nodes(),x,y)));
				return 1;									//returns 1 if an M augmenting path is found
			}
			
			else if(T[y]==0){
				Tree.add_edge(x,y);
				T[(M.node(y)).neighbors().at(0)]=1;
				T_even.push_back((M.node(y)).neighbors().at(0));
				T[y]=2;
				return 0;									//returns 0 if an edge was added to T
			}
			
			else {return 0;}//TODO: shrink and update		//returns 0 if an odd cycle was found
		}	
	}
	return 2;												//returns 2 if the tree is frustrated
}

int build_up_tree(ED::Graph const & G, ED::Graph & M, ED::NodeId start) {
	ED::Graph Tree{G.num_nodes()};
	std::vector <int> T(G.num_nodes());
	std::vector <ED::NodeId> T_even;
	
	//to simulate the shrinking we use partition classes
	std::vector <int> Partition(G.num_nodes());
	std::vector <std::vector <ED::NodeId>> Partition_classes(G.num_nodes());
	
	for(unsigned int i=0; i<G.num_nodes(); i++){
		Partition.at(i)=i;		//at the beginning each vertex has its own partition class 
		Partition_classes.at(i).push_back(i);
	}
	
	T_even.push_back(start);
	T[start]=1;
	
	int i=0;

	while(i==0){
		i=grow_tree(G,M,start, Tree, T, T_even, Partition, Partition_classes);
	}
	if(i=1){
		if(check_M_perfect(M)==1){
			return 1;
		}
		else return 0;
	}
	else return 2;
}

ED::Graph perfect_Matching(ED::Graph const &G) {
	ED::Graph M{G.num_nodes()};
	ED::NodeId start=0;
	if(G.num_nodes()==0 || G.num_nodes()==1){ 
		std::cout<< "G is empty or has only one vertex";
		return M;
	}
	int i=build_up_tree(G, M, start);
	if(i==1) {
		std::cout<< "M has perfect matching";
		return M;
	}
	else if(i==2) {
		std::cout<<"M has no perfect matching";
		return M;
	}
	
	
	while(i==0) {
		find_M_exposed_node(M, start);
		i=build_up_tree(G,M,start);
	}
	
	if(i==1) {
		std::cout<< "M has perfect matching";
		return M;
	}
	else if(i==2) {
		std::cout<<"M has no perfect matching";
		return M;
	}
	
	return M;
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
