#include <iostream>
#include <cstdlib>
#include <utility>

//#include "graph.cpp"  //already contained in "readGraph"
//#include "graph.hpp"
#include "readGraph.cpp"
#include "UnionFind.cpp"


void agument_M(ED::Graph & M, std::vector <ED::NodeId> const & path) {
	for (unsigned int j=1; j < path.size()-1; j++) {
		M.delete_edge(path[j], path[j+1]);
		j++;
	}
	for (unsigned int j=0; j < path.size()-1; j++) {
		M.add_edge(path[j],path[j+1]);
		j++;
	}
}

//Only for testing, ToDo: Delete
//checks if "path" is M-augmenting
bool path_M_augmenting(ED::Graph const & M, std::vector <ED::NodeId> const & path) {
	int i=0;
	for(unsigned int j=0; j< path.size(); j++){
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

//Warning: runtime OK but dont work correctly
bool find_M_agumenting_path_inclass(ED::Graph const & Tree, ED::Graph const & M, std::vector <bool> & visited, UF::UnionFind & Partition, std::vector <ED::NodeId> & p, ED::NodeId start, ED::NodeId end) {
	
	visited.at(start)=1;
	p.push_back(start);
	
	for(auto const & neighbor_id : (Tree.node(start)).neighbors()) {
		if(M.node(neighbor_id).neighbors().at(0)==end){
			p.push_back(neighbor_id);
			p.push_back(M.node(neighbor_id).neighbors().at(0));
			return 1;
		}
		else if(Partition.equals(neighbor_id, start) && visited.at(neighbor_id)!=0) {
			if(Partition.equals(M.node(neighbor_id).neighbors().at(0),start)) {
				p.push_back(neighbor_id);
				if(find_M_agumenting_path_inclass(Tree,M,visited,Partition,p,M.node(neighbor_id).neighbors().at(0),end)==1) return 1;
				p.pop_back();
				p.pop_back();
			}
		}
	}
	return 0;
}

//Warning: dont work correctly
std::vector <ED::NodeId> find_M_agumenting_path(ED::Graph const & M, ED::Graph const & Tree, UF::UnionFind & Partition, std::vector <ED::NodeId> const & p) {
	int start=0;
	int end=0;
	int k=Partition.find(p.at(0));
	std::vector <bool> visited(M.num_nodes(), 0);
	
	std::vector <ED::NodeId> path;
	path.push_back(p.at(0));
	
	for(unsigned int j=1; j<p.size(); j++) {
		if(k==Partition.find(p.at(j))){
			end=j;
		}
		else {
			if(end!=start){
				find_M_agumenting_path_inclass(Tree, M, visited, Partition, path, p.at(start), p.at(end));
			}
			else path.push_back(p.at(j));
			k=Partition.find(p.at(j));
			start=j;
			end=j;
		}
	}
	return path;
}

bool DFSI(ED::Graph const & Tree, ED::NodeId v, ED::NodeId w, std::vector <ED::NodeId> & path, std::vector <bool> & visited, bool i){
	visited[v]=1;
	auto const & node = Tree.node(v);
	path.push_back(v);
	
	for (auto const & neighbor_id : node.neighbors()) {
		if(neighbor_id==w) {
			path.push_back(w);
			return 1;
		}
		else if(visited[neighbor_id]==0){
			if(DFSI(Tree, neighbor_id, w, path, visited,1-i)) return 1;
			path.pop_back();
		}
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

//checks whether v is covered by the matching M
bool check_node_M_exposed(ED::NodeId v, ED::Graph const & M) {
	auto const & node = M.node(v);
	if (node.neighbors().size()!=0) {
		return 0;
	}
	return 1;
}

ED::NodeId find_M_exposed_node(ED::Graph const & M){
	for(unsigned int i=0; i<M.num_nodes(); i++){
		if((M.node(i)).neighbors().size()==0){
			return i;
		}
	}
	return M.num_nodes();
}

//checks whether M is a perfect matching
bool check_M_perfect(ED::Graph const & M){
	for (int node_id = 0; node_id < (int)M.num_nodes(); ++node_id) {
		if (((M.node(node_id)).neighbors()).size()==0) {
			return 0;
		}
	}
	return 1;
}

void update_partition_edges(ED::Graph G, std::vector <ED::NodeId> path, UF::UnionFind & Partition,std::vector <std::pair<ED::NodeId, ED::NodeId>> Edges){
	
	//ToDo update k at "Partition.unite(k,path.at(i))"
	
	int k=Partition.find(path.at(0));
	int current=k;
	
	//Position in Tree, 0 for even 1 for odd
	int position=0;
	
	if(path.size()>1) {
		for(unsigned int i=1; i<path.size(); i++){
			//if the next vertex on "path" lies in another partitionclass, we unite the current class and class k  
			if(Partition.find(path.at(i))!=current){
				//position is 0 if the previous class was even, so the current vertex is odd
				if(position==0) {
					//update "Edges"
					auto const & node = G.node(i);
					for (auto const & neighbor_id : node.neighbors()) {
						Edges.push_back(std::make_pair(i,neighbor_id));
					}
				}
				position=1-position;
				//update "Partition"
				Partition.unite(k,path.at(i));
				current=Partition.find(path.at(i));
			}
		}
	}
	//change the position of the partition in tree to even
	Partition.set_position_even(k);
}


//grow tree corresponds to one iteration of the while loop in the perfect matching algorithm
int grow_tree (ED::Graph const & G, ED::Graph & M, ED::NodeId const root, ED::Graph & Tree, UF::UnionFind & Partition, ED::NodeId const start, ED::NodeId end, std::vector <std::pair<ED::NodeId, ED::NodeId>> & Edges) 
{	ED::NodeId even;
	if(Partition.get_position_in_tree(end)==0 && check_node_M_exposed(end,M)){
		Tree.add_edge(end,start);
		Partition.set_position_odd(start);
		//agument_M(M,find_M_agumenting_path(M,Tree,Partition,DFS(Tree, Tree.num_nodes(),root,end))); // ToDO:Implement function correctly
		agument_M(M,DFS(Tree, Tree.num_nodes(),root,end));
		return 1;													//returns 1 if an M augmenting path is found
	}
			
	else if(Partition.get_position_in_tree(end)==0){
		//update Tree and Positions in Tree
		even=(M.node(end)).neighbors().at(0);
		Tree.add_edge(start,end);  
		Tree.add_edge(end,even);
		Partition.set_position_even(even);
		Partition.set_position_odd(end);
				
		//update "Edges" by adding all edges incident to the vertex "even" 
		auto const & node = G.node(even);
		for (auto const & neighbor_id : node.neighbors()) {
			if(neighbor_id!=end) Edges.push_back(std::make_pair(even,neighbor_id));
		}
		return 0;													//returns 0 if an edge was added to T
	}
			
	else {
		//update "Partition" and update "Edges" by adding every edge incident to an odd vertex 
		//Warning: we do not check if this edge already exists
		//to avoid double edges we could check if start and end have the same Partitionclass
		update_partition_edges(G,DFS(Tree, Tree.num_nodes(),start,end), Partition, Edges);
		Tree.add_edge(start,end);
		return 0;													//returns 0 if an odd cycle was found
	}		
	return 2;														//returns 2 if the tree is frustrated
}



//builds up a tree and returns 0 if M is not perfect, 1 if M is alreday perfect and 2 if the tree is frustrated
int build_up_tree(ED::Graph const & G, ED::Graph & M) {
	ED::Graph Tree{G.num_nodes()};
	UF::UnionFind Partition(G.num_nodes());
	
	//"Edges" is the set of edges incident to at least one vertex in the tree
	//Warning: "Edges" may contain the same edge two times
	std::vector <std::pair<ED::NodeId, ED::NodeId>> Edges;
	
	ED::NodeId root,start,end;
	int i=0;

	root=find_M_exposed_node(M);
	
	if(root==M.num_nodes()) {
		return 1;
	}
	
	//build up tree with this root
	else {
		auto const & node = G.node(root);
		for (auto const & neighbor_id : node.neighbors()){
			Edges.push_back(std::make_pair(root, neighbor_id));
		}
		
		while(Edges.size()>0 && i==0) {
			start=Edges.at(Edges.size()-1).first;
			end=Edges.at(Edges.size()-1).second;
			//if "end" has odd distance to the root, we don't consider this edge
			//if "end" becomes even after some shrinking operation, this edge will be added to "Edges" during the shrinking
			if(Partition.get_position_in_tree(end)!=2) i=grow_tree(G,M,root,Tree,Partition,start,end,Edges); 
		}
		
		if(i==1){
			if(check_M_perfect(M)) return 1;
			else return 0;
		}
	}
	return 2;
}


ED::Graph perfect_Matching(ED::Graph const &G) {
	ED::Graph M{G.num_nodes()};
	if(G.num_nodes()==0 || G.num_nodes()==1){ 
		std::cout<< "G is empty or has only one vertex";
		return M;
	}
	int i=0;

	
	while(i==0) {
		i=build_up_tree(G,M);
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

//testing
int main()
{	
	ED::Graph g1{8};
	ED::Graph g2{6};
	ED::Graph g3{10};
	ED::Graph g4{10};
	
	g1.add_edge(6,0);
	g1.add_edge(0,1);
	g1.add_edge(2,3);
	g1.add_edge(4,5);
	g1.add_edge(1,2);
	g1.add_edge(2,7);
	g1.add_edge(4,3);
	g1.add_edge(3,5);
	g1.add_edge(1,4);
	
	g2.add_edge(1,2);
	g2.add_edge(3,0);
	g2.add_edge(4,5);
	g2.add_edge(0,1);
	g2.add_edge(2,3);
	g2.add_edge(0,4);
	
	g3.add_edge(0,9);
	g3.add_edge(0,1);
	g3.add_edge(1,2);
	g3.add_edge(2,3);
	g3.add_edge(3,4);
	g3.add_edge(4,5);
	g3.add_edge(5,6);
	g3.add_edge(6,7);
	g3.add_edge(7,8);
	g3.add_edge(8,1);
	

	std::cout<< perfect_Matching(g1);
	std::cout<< perfect_Matching(g2);
	std::cout<< perfect_Matching(g3);

   return 1;
}
