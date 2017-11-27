#include "union_find.hpp"

#include<vector>
#include<stdexcept>
#include<cstddef>

namespace UF
{
/*
Initialize the union find structure
the vector forest represents the trees, where forest[i] = j means that vertex i has j as a parent, and all roots point to themselves
Upon initalization, each number is it's own partition class, meaning each number is a root
All ranks are zero
*/
UnionFind::UnionFind(std::size_t size)
	: _forest(size),
	_ranks(size, 0),
	_size(size)
	{
		for (unsigned int i = 0; i<size; i++){
			_forest.at(i).push_back(i);
			//Position in tree
			_forest.at(i).push_back(0);
			//apex of pseudonode
			//warning: position and apex only updated for the root of a pseudonode
			_forest.at(i).push_back(-1);
		}
	}


/*
Find the root of the partition class to which number i belongs
Checks if index i is in {0, ..., n-1}, throws std::out_of_range otherwise
Returns the root of the partition class
*/
int UnionFind::find(unsigned int i){
	if (i >= _size){
		throw std::out_of_range("UnionFind::find");
	}
	int current = i;
	int parent = _forest.at(i).at(0);
	while (parent != current){
		current = parent;
		parent = _forest.at(current).at(0);
	}
	_forest.at(i).at(0)=current;
	return current;
}


int UnionFind::get_position_in_tree(int i) {
	return _forest.at(find(i)).at(1);
}

void UnionFind::set_position_even(int i) {
	_forest.at(find(i)).at(1)=1;
}
void UnionFind::set_position_odd(int i) {
	_forest.at(find(i)).at(1)=2;
}
/*
Checks if i and j belong to the same partition class or not
*/
bool UnionFind::equals(int i, int j){
	return (find(i)==find(j));
}

/*
Unites the partition classes to which i and j belong
If i and j belong to the same class, a std::logic_error is thrown
Internally, the smaller of the trees i and j belong to is attached to the root of the bigger one
If neccessary, the rank of the roots of the trees are incremented
*/
void UnionFind::unite(int i, int j){
	int root_i = find(i);
	int root_j = find(j);
	if (root_i == root_j){
		throw std::logic_error("UnionFind::unite cannot be executed for two equal partition classes.");
	}
	if (_ranks[root_i] > _ranks[root_j]){
		_forest.at(root_j).at(0) = root_i;
	}
	else {
		_forest.at(root_i).at(0) = root_j;
		if (_ranks[root_i] == _ranks[root_j]){
			_ranks[root_j] += 1;
		}
	}
}
}//namespace UF

