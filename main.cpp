#include <iostream>
#include <cstdlib>
#include <utility>

#include <tuple> //std::tie
#include "read_graph.hpp"
#include "union_find.hpp"

//augment M along path
void augment_M(ED::Graph & M, std::vector <ED::NodeId> const & path) {
	for (unsigned int j=1; j < path.size()-1; j++) {
		M.delete_edge(path[j], path[j+1]);
		j++;
	}
	for (unsigned int j=0; j < path.size()-1; j++) {
		M.add_edge(path[j],path[j+1]);
		j++;
	}
}
//checks whether v is covered by the matching M
bool check_node_M_exposed(ED::NodeId v, ED::Graph const & M) {
	auto const & node = M.node(v);
	if (node.neighbors().size()!=0) {
		return 0;
	}
	return 1;
}
//checks whether v is M-exposed
ED::NodeId find_M_exposed_node(ED::Graph const & M, std::vector<bool> const & considered){
	for(unsigned int i=0; i<M.num_nodes(); i++){
		if(considered.at(i) && (M.node(i)).neighbors().size()==0){
			return i;
		}
	}
	return M.num_nodes();
}

//checks whether M is a perfect matching
bool check_M_perfect(ED::Graph const & M){
	for (unsigned int node_id = 0; node_id <M.num_nodes(); ++node_id) {
		if (((M.node(node_id)).neighbors()).size()==0) {
			return 0;
		}
	}
	return 1;
}

// DFS only computes a path in T, it has not to be M augmenting in the partition classes
// so we compute in each partition class the M augmenting path

bool find_M_augmenting_path_inclass(ED::Graph const & Tree, ED::Graph const & M, std::vector <bool> & visited, UF::UnionFind & Partition, std::vector <ED::NodeId> & p, ED::NodeId start, ED::NodeId end) {
	//visit start
	visited.at(start)=1;
	//add start to path
	p.push_back(start);
	
	//consider the neighbors of start
	for(auto const & neighbor_id : (Tree.node(start)).neighbors()) {
		
		if(!check_node_M_exposed(neighbor_id,M)&& neighbor_id!=end){
			//always only consider the vertex mached to neighbor_id to make sure that the path is M-augmenting
			
			if(M.node(neighbor_id).neighbors().at(0)==end){
				//M-augmenting path was found
				p.push_back(neighbor_id);
				p.push_back(M.node(neighbor_id).neighbors().at(0));
				return 1;
			}
			
			//make sure that neighbor_id don't leave the partitionclass and the vertex mached to it was not visited so far 
			else if(Partition.equals(neighbor_id, start) && visited.at(M.node(neighbor_id).neighbors().at(0))==0) {
				
				//make sure that the mached vertex don't leave the partitionclass
				if(Partition.equals(M.node(neighbor_id).neighbors().at(0),start)) {
					p.push_back(neighbor_id);
					//neighbor_id is marked as visited to avoid circuits in the path
					visited.at(neighbor_id)=1;
					//call this function again with new start: M.node(neighbor_id).neighbors().at(0)
					if(find_M_augmenting_path_inclass(Tree,M,visited,Partition,p,M.node(neighbor_id).neighbors().at(0),end)==1) return 1;
					
					//if end wasn't found
					//mark neighbor_id as not visited
					visited.at(neighbor_id)=0;
					p.pop_back();
					p.pop_back();
				}
			}
		}
	}
	return 0;
}
//Modifies the path p sucht that p is M-augmenting in each pseudonode whith non trivial intersection to p
std::vector <ED::NodeId> find_M_augmenting_path(ED::Graph const & M, ED::Graph const & Tree, UF::UnionFind & Partition, std::vector <ED::NodeId> const & p) {
	//p(start) is the first vertex of the current pseudonode lying on p
	int start=0;
	//p(end) is the last vertex of the current pseudonode lying on p
	int end=0;
	//current pseudonode
	int current=Partition.find(p.at(0));
	
	std::vector <bool> visited(M.num_nodes(), 0);
	std::vector <ED::NodeId> path;
	
	for(unsigned int j=1; j<p.size(); j++) {
		
		if(current==Partition.find(p.at(j))){
			//update end
			end=j;
		}
		else {
			
			if(end!=start){ //non trivial intersection of p and pseudonode
				
				if(!check_node_M_exposed(p.at(end),M)) {
					
					//find an M-augmenting path in the pseudonode current from p(start) to p(end)
					if(Partition.equals(M.node(p.at(end)).neighbors().at(0),p.at(end)))
						find_M_augmenting_path_inclass(Tree, M, visited, Partition, path, p.at(start), p.at(end));
					//find an M-augmenting path in the pseudonode current from p(end) to p(start)
					else 
						find_M_augmenting_path_inclass(Tree, M, visited, Partition, path, p.at(end), p.at(start));
				}
				
				//find an M-augmenting path in the pseudonode current from p(end) to p(start)
				else 
					find_M_augmenting_path_inclass(Tree, M, visited, Partition, path, p.at(end), p.at(start));
			}
			
			else path.push_back(p.at(j-1));	//trivial intersection of p and pseudonode
			//update of current, start and end
			current=Partition.find(p.at(j));
			start=j;
			end=j;
		}
	}
	
	//add last vertex of p to path
	path.push_back(p.at(p.size()-1));
	return path;
}

bool DFSI(ED::Graph const & Tree, ED::NodeId v, ED::NodeId w, std::vector <ED::NodeId> & path, std::vector <bool> & visited, bool i){
	//visit v
	visited[v]=1;
	auto const & node = Tree.node(v);
	//add v to path
	path.push_back(v);
	
	//consider all neighbors of v
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


void update_partition_edges(ED::Graph const & G, std::vector<bool> const & considered, std::vector <ED::NodeId> const & path, UF::UnionFind & Partition, std::vector <std::pair<ED::NodeId, ED::NodeId>> & Edges){
	
	int k=Partition.find(path.at(0));
	int current=k;
	
	if(path.size()>1) {
		for(unsigned int i=1; i<path.size(); i++){
			//if the next vertex on "path" lies in another partitionclass, we unite the current class and class k  
			if(Partition.find(path.at(i))!=current){
				if(Partition.get_position_in_tree(path.at(i))==2) {
					//update "Edges"
					auto const & node = G.node(path.at(i));
					for (auto const & neighbor_id : node.neighbors()) {
						if (considered.at(neighbor_id)&& neighbor_id!=path.at(i-1)&& neighbor_id!=path.at(i+1)){
							Edges.push_back(std::make_pair(path.at(i),neighbor_id));
						}
					}
				}
				//update "Partition"
				Partition.unite(k,path.at(i));
				current=Partition.find(path.at(i));
			}
		}
	}
	//change the position of the partition in tree to even
	Partition.set_position_even(k);
}


//grow tree corresponds to one iteration of the while loop in the perfect matching algorithm, the considered edge  is {start,end}
int grow_tree (ED::Graph const & G, ED::Graph & M, std::vector<bool> const & considered, ED::NodeId const root, ED::Graph & Tree, UF::UnionFind & Partition, ED::NodeId const start, ED::NodeId end, std::vector <std::pair<ED::NodeId, ED::NodeId>> & Edges) 
{	ED::NodeId even;
	//if end is  M-exposed, augment M along a path from root to end
	if(Partition.get_position_in_tree(end)==0 && check_node_M_exposed(end,M)){
		Tree.add_edge(end,start);
		Partition.set_position_odd(end);
		//first: find a path from root to end in Tree
		//second: make it M-augmenting and then augment M
		augment_M(M,find_M_augmenting_path(M,Tree,Partition,DFS(Tree, Tree.num_nodes(),root,end)));
		return 1;		//returns 1 if M-augmenting path was found
	}
	
	//if end is not part of Tree but not M-exposed			
	else if(Partition.get_position_in_tree(end)==0){
		//update Tree and Positions in Tree
		even=(M.node(end)).neighbors().at(0);
		if (considered.at(even)){
			Tree.add_edge(start,end);  
			Tree.add_edge(end,even);
			Partition.set_position_even(even);
			Partition.set_position_odd(end);
				
			//update "Edges" by adding all edges incident to the vertex "even" 
			auto const & node = G.node(even);
			for (auto const & neighbor_id : node.neighbors()) {
				if(neighbor_id!=end && considered.at(neighbor_id)) Edges.push_back(std::make_pair(even,neighbor_id));
			}
		}
		return 0;													//returns 0 if an edge was added to T
	}
			
	//if odd cycle was found	
	else if(!Partition.equals(start,end)) {
		//update "Partition" and update "Edges" by adding every edge incident to an odd vertex 
		//to avoid double edges we could check if start and end have the same Partitionclass
		update_partition_edges(G, considered, DFS(Tree, Tree.num_nodes(),start,end), Partition, Edges);
		Tree.add_edge(start,end);
		return 0;													//returns 0 if an odd cycle was found
	}		
	return 0;														
}



//builds up a tree and returns 0 if M is not maximum but the tree is not frustrated, 1 if M is alreday maximum and 2 if the tree is frustrated
int build_up_tree(ED::Graph const & G, ED::Graph & M, std::vector<bool> & considered) {
	ED::Graph Tree{G.num_nodes()};
	UF::UnionFind Partition(G.num_nodes());
	
	//"Edges" is the set of edges incident to at least one vertex in the tree
	//Warning: "Edges" may contain the same edge two times
	std::vector <std::pair<ED::NodeId, ED::NodeId>> Edges;
	
	ED::NodeId root,start,end;
	int i=0;

	root=find_M_exposed_node(M, considered);
	if(root==M.num_nodes()) {
		return 1;
	}
	
	//build up tree with this root
	else {
		//update Position for root
		Partition.set_position_even(root);
		auto const & node = G.node(root);
		for (auto const & neighbor_id : node.neighbors()){
			if (considered.at(neighbor_id)) {
				//add all edges incident to "root"
				Edges.push_back(std::make_pair(root, neighbor_id));
			}
		}
		
		while(Edges.size()>0 && i==0) {
			start=Edges.at(Edges.size()-1).first;
			end=Edges.at(Edges.size()-1).second;
			Edges.pop_back();
			//if "end" has odd distance to the root, we don't consider this edge
			//if "end" becomes even after some shrinking operation, this edge will be added to "Edges" during the shrinking
			if(Partition.get_position_in_tree(end)!=2) i=grow_tree(G,M,considered, root,Tree,Partition,start,end,Edges); 
		}
		
		if(i==1){
			if(check_M_perfect(M)) return 1;
			else return 0;
		}
	}
	//"remove" the frustrated tree by setting it's vertices to not considered
	for (unsigned int j=0; j<Tree.num_nodes(); j++){
		auto const & node = Tree.node(j);
		if (node.degree() != 0||j==root){
			considered.at(j) = false;
		}
	}
	return 2;
}

/*
Construct a maximum matching in G using the hint matching M
*/
ED::Graph maximum_matching(ED::Graph const &G, ED::Graph &M) {
	//in the first iteration, all vertices are considered
	std::vector<bool> considered(G.num_nodes(),true);

	if(G.num_nodes()==0 || G.num_nodes()==1){ 
		return ED::Graph{G.num_nodes()};
	}
	int i=0;

	
	while(i==0|| i==2) {
		i=build_up_tree(G, M, considered);
	}
	
	if(i==1) {
		return M;
	}
	
	return M;
}


int main(int argc, char* argv[])
{	
	if (argc == 2){
		int success;
		ED::Graph graph(0);
		std::tie(success, graph) = ED::read_graph(argv[1]);
		if (success != EXIT_SUCCESS){
			std::cerr << "The graph could not be read successfully.";
			return EXIT_FAILURE;
		}
		ED::Graph no_hint_matching(graph.num_nodes());
		ED::Graph matching = maximum_matching(graph, no_hint_matching);
		std::cout << "The following maximum matching was found:\n";
		std::cout << matching;
		return EXIT_SUCCESS;
	}
	else if (argc == 3){
		int success_1, success_2;
		ED::Graph graph(0);
		ED::Graph hint_matching(0);
		std::tie(success_1, graph) = ED::read_graph(argv[1]);
		std::tie(success_2, hint_matching) = ED::read_graph(argv[2]);
		if (success_1 != EXIT_SUCCESS || success_2 != EXIT_SUCCESS){
			std::cerr << "The graph or hint matching could not be read successfully.";
			return EXIT_FAILURE;
		}
		ED::Graph matching = maximum_matching(graph, hint_matching);
		std::cout << "The following maximum matching was found:\n";
		std::cout << matching;
		return EXIT_SUCCESS;
	}
	else{
		std::cout << "Wrong number of arguments. One argument filename (file enconding the graph) is obligatory, second argument for hint matching is optional, no futher arguments can be processed.";
		return EXIT_FAILURE;
	}
}
