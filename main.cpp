#include <iostream>
#include <cstdlib>
#include <utility>

#include "readGraph.cpp"
#include "UnionFind.cpp"



//ugment M along path
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
ED::NodeId find_M_exposed_node(ED::Graph const & M){
	for(unsigned int i=0; i<M.num_nodes(); i++){
		if((M.node(i)).neighbors().size()==0){
			return i;
		}
	}
	return M.num_nodes();
}

//ToDo remove, only for testing
bool check_path_simple(std::vector <ED::NodeId>const &path, unsigned int size) {
	std::vector <bool> visited(size, 0);
	for(unsigned int i=0; i<path.size(); i++) {
		if(visited.at(path.at(i))==0) {visited.at(path.at(i))=1;}
		else return 0;
	}
	return 1;
}

//checks whether M is a perfect matching
bool check_M_perfect(ED::Graph const & M){
	for (int node_id = 0; node_id < (int)M.num_nodes(); ++node_id) {
		if (((M.node(node_id)).neighbors()).size()!=1) {
			return 0;
		}
	}
	return 1;
}

/* DFS only computes a path in T, it has not to be M augmenting in the partition classes
 so we compute in each partition class the M augmenting path*/
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


bool DFSI(ED::Graph const & Tree, ED::NodeId v, ED::NodeId w, std::vector <ED::NodeId> & path, std::vector <bool> & visited){
	//visit v
	visited[v]=1;
	auto const & node = Tree.node(v);
	
	//add v to path
	path.push_back(v);
	
	//consider all neighbors of v
	for (auto const & neighbor_id : node.neighbors()) {
		
		if(neighbor_id==w) {
			//path was found
			path.push_back(w);
			return 1;
		}
		
		else if(visited[neighbor_id]==0){
			//call function with new v: neighbor_id
			if(DFSI(Tree, neighbor_id, w, path, visited)) return 1;
			//if path from beighbor_id to w wasn't found
			path.pop_back();
		}
	} 
	return 0;
}

//Depth first search
std::vector <ED::NodeId> DFS(ED::Graph const & T, const int size,  ED::NodeId v, ED::NodeId w) {
	std::vector <ED::NodeId> path;
	std::vector <bool> visited(size, 0);
	//find path from v to w
	DFSI(T,v,w,path,visited);
	return path;
}



//update "Partition" and "Edges"
void update_partition_edges(ED::Graph const &G, std::vector <ED::NodeId> const & path, UF::UnionFind & Partition, std::vector <std::pair<ED::NodeId, ED::NodeId>> & Edges){
	
	int k=Partition.find(path.at(0));
	int current=k;
	
	if(path.size()>1) {
		for(unsigned int i=1; i<path.size(); i++){
			//if the next vertex on "path" lies in another partitionclass, we unite the current class and class k  
			if(Partition.find(path.at(i))!=current){
				//position is 0 if the previous class was even, so the current vertex is odd
				if(Partition.get_position_in_tree(path.at(i))==2) {
					//update "Edges"
					auto const & node = G.node(path.at(i));
					for (auto const & neighbor_id : node.neighbors()) {
						if(neighbor_id!=path.at(i-1)&& neighbor_id!=path.at(i+1)){
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
int grow_tree (ED::Graph const & G, ED::Graph & M, ED::NodeId const root, ED::Graph & Tree, UF::UnionFind & Partition, ED::NodeId const start, ED::NodeId const end, std::vector <std::pair<ED::NodeId, ED::NodeId>> & Edges) 
{	ED::NodeId even;
	//if end is  M-exposed, augment M along a path from root to end
	if(check_node_M_exposed(end,M) && Partition.get_position_in_tree(end)==0){
		Tree.add_edge(end,start);
		Partition.set_position_odd(end);
		//first: find a path from root to end in Tree
		//second: make it M-augmenting and then augment M
		//if(check_path_simple(find_M_augmenting_path(M,Tree,Partition,DFS(Tree, Tree.num_nodes(),root,end)),M.num_nodes())){
			augment_M(M,find_M_augmenting_path(M,Tree,Partition,DFS(Tree, Tree.num_nodes(),root,end)));
			return 1;
		//}
		//agument_M(M,DFS(Tree, Tree.num_nodes(),root,end));
		//return 2;													//returns 1 if an M augmenting path is found
	}
	
	//if end is not part of Tree but not M-exposed	
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
	
	//if odd cycle was found	
	else if(!Partition.equals(start,end)){
		//update "Partition" and update "Edges" by adding every edge incident to an odd vertex 
		//to avoid double edges we could check if start and end have the same Partitionclass
		update_partition_edges(G,DFS(Tree, Tree.num_nodes(),start,end), Partition, Edges);
		Tree.add_edge(start,end);
		return 0;													//returns 0 if an odd cycle was found
	}		
	return 0;
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
		//update Position for root
		Partition.set_position_even(root);
		
		auto const & node = G.node(root);
		for (auto const & neighbor_id : node.neighbors()){
			//add all edges incident to "root"
			Edges.push_back(std::make_pair(root, neighbor_id));
		}
		
		while(Edges.size()>0 && i==0) {
			start=Edges.at(Edges.size()-1).first;
			end=Edges.at(Edges.size()-1).second;
			Edges.pop_back();
			
			//if "end" has odd distance to the root, we don't consider this edge
			//if "end" becomes even after some shrinking operation, this edge will be added to "Edges" during the shrinking
			if(Partition.get_position_in_tree(end)!=2) i=grow_tree(G,M,root,Tree,Partition,start,end,Edges); 
		}
		
		if(i==1){
			if(check_M_perfect(M)) return 1;
			else return 0;
		}
	}
	//ToDo remove tree
	return 2;
}
//ToDo delete, only for testing
bool check_M_part_of_G(ED::Graph M, ED::Graph G) {
	bool check;
	for(unsigned int i=0; i<M.num_nodes();i++) {
		for(unsigned int k=0;k<G.node(i).neighbors().size();k++) {
			if(M.node(i).neighbors().at(0)==G.node(i).neighbors().at(k)) check=1;
		}
		if(check==0) return 0;
	}
	return 1;
}

//perfect Matching algorithm
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
{	std::cout<<perfect_Matching(ED::read_graph("lu980.dmx").second);
	//std::cout<<check_M_perfect(perfect_Matching(ED::read_graph("ei8246.dmx").second))<<"\n";
	//std::cout<<check_M_part_of_G(perfect_Matching(ED::read_graph("ei8246.dmx").second),ED::read_graph("ei8246.dmx").second);
	
	ED::Graph g{8};
	ED::Graph g1{8};
	g.add_edge(0,3);
	g.add_edge(2,1);
	g.add_edge(2,0);
	g.add_edge(0,1);
	
	
	g1.add_edge(6,0);g1.add_edge(1,4);
	g1.add_edge(1,2);
	
	g1.add_edge(5,3);
	g1.add_edge(2,7);
	
	g1.add_edge(0,1);
	g1.add_edge(2,3);
	g1.add_edge(4,5);
	std::cout<<perfect_Matching(g1);std::cout<<check_M_part_of_G(g,g1);
   return 1;
}
