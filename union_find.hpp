#ifndef UNION_FIND_HPP
#define UNION_FIND_HPP

#include <cstddef> // std::size_t
#include <vector>


namespace UF{

/*
Class for UnionFind structure
on the numbers {0, ..., n-1}
in order to store a partition of this set
*/
class UnionFind
{
public:
	//On initialization, each number forms it's own partition class
	UnionFind(std::size_t size);

	//Find the representative (root) of i's partition class
	int find(unsigned int i);

	//Compare if i and j belong to the same partition class
	bool equals(int i, int j);

	//Unite the partition classes of i and j
	void unite(int i, int j);
	
	int get_position_in_tree(int i);
	
	void set_position_even(int i);
	
	void set_position_odd(int i);

private:
	//Internally save the partition classes as trees in forest
	std::vector<std::vector <int>> _forest;

	//Height of every vertex in the respective tree
	//@warning: Note that ranks are only kept up to date for roots of trees!
	std::vector<int> _ranks;

	//size of the set {0, ..., n-1}
	std::size_t _size;
};

} // namespace UF

#endif
