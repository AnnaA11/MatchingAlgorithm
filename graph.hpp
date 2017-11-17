#ifndef GRAPH_HPP
#define GRAPH_HPP

/**
   @file graph.hpp

   @brief This file provides a simple class @c Graph to model unweighted undirected graphs.
**/

#include <cstddef> // std::size_t
#include <iosfwd> // std::ostream fwd declare
#include <limits>
#include <vector>

namespace ED // for Edmonds
{

using size_type = std::size_t;

//! Always use these typedefs to identify nodes by their numbers.
//! Use different typedefs for internal standard indexing starting with 0
//! and DIMACS-based indexing starting with 1.
//! Even better, one could use strong typedefs (i.e. use an enum class
//! or a custom struct providing cast operator to \c size_type s.t.
//! mixup of index types is avoided.
//! This is not done here for simplicity.
using NodeId = size_type;
using DimacsId = size_type;

/** Useful constant different from the id of any actual node: **/
NodeId constexpr invalid_node_id = std::numeric_limits<NodeId>::max();
DimacsId constexpr invalid_dimacs_id = std::numeric_limits<DimacsId>::max();

/**
   Nodes in DIMACS files are counted from 1, but here we count them from 0 so they match their std::vector indices.
   These two trivial functions should help make the transition between the two models clear (instead of just having
   some unexplained -1's and +1's in the middle of the code.
**/
NodeId from_dimacs_id(DimacsId const dimacs_id); //!< Subtracts 1 (throws if @c dimacs_id is 0)
DimacsId to_dimacs_id(NodeId const node_id);     //!< Adds 1 (throws if overflow would occur)

/**
   @class Node

   @brief A @c Node stores an array of neighbors (via their ids).

   @note The neighbors are not necessarily ordered, so searching for a specific neighbor takes O(degree)-time.
**/
class Node
{
public:
   typedef std::size_t size_type;

   /** @brief Create an isolated node (you can add neighbors later). **/
   Node() = default;

   /** @return The number of neighbors of this node. **/
   size_type degree() const;

   /** @return The array of ids of the neighbors of this node. **/
   std::vector<NodeId> const & neighbors() const;

private:
   friend class Graph;

   /**
      @brief Adds @c id to the list of neighbors of this node.
      @warning Does not check whether @c id is already in the list of neighbors (a repeated neighbor is legal, and
      models parallel edges).
      @warning Does not check whether @c id is the identity of the node itself (which would create a loop!).
   **/
   void add_neighbor(NodeId const id);
   void delete_neighbor();

   std::vector<NodeId> _neighbors;
}; // class Node

/**
   @class Graph

   @brief A @c Graph stores an array of @c Node s, but no array of edges. The list of edges is implicitly given
   by the fact that the nodes know their neighbors.

   This class models undirected graphs only (in the sense that the method @c add_edge(node1, node2) adds both @c node1
   as a neighbor of @c node2 and @c node2 as a neighbor of @c node1). It also forbids loops, but parallel edges are
   legal.

   @warning Nodes are numbered starting at 0, as is usually done in programming,
    instead starting at 1, as is done in the DIMACS format that your program should take as input!
    Be careful.
**/
class Graph
{
public:
   typedef std::size_t size_type;

   /**
      @brief Creates a @c Graph with @c num_nodes isolated nodes.

      The number of nodes in the graph currently cannot be changed. You can only add edges between the existing nodes.
   **/
   Graph(NodeId const num_nodes);

   /** @return The number of nodes in the graph. **/
   NodeId num_nodes() const;

   /** @return The number of edges in the graph. **/
   size_type num_edges() const;

   /**
      @return A reference to the id-th entry in the array of @c Node s of this graph.
   **/
   Node const & node(NodeId const id) const;

   /**
      @brief Adds the edge <tt> {node1_id, node2_id} </tt> to this graph.

      Checks that @c node1_id and @c node2_id are distinct and throws an exception otherwise.
      This method adds both @c node1_id as a neighbor of @c node2_id and @c node2_id as a neighbor of @c node1_id.

      @warning Does not check that the edge does not already exist, so this class can be used to model non-simple graphs.
   **/
   void add_edge(NodeId node1_id, NodeId node2_id);
   void delete_edge(NodeId a, NodeId b);

   /**
     @brief Prints the graph to the given ostream in DIMACS format.
   **/
   friend std::ostream & operator<<(std::ostream & str, Graph const & graph);

private:
   std::vector<Node> _nodes;
   std::size_t _num_edges;
}; // class Graph

//BEGIN: Inline section

inline
Node::size_type Node::degree() const
{
   return neighbors().size();
}

inline
std::vector<NodeId> const & Node::neighbors() const
{
   return _neighbors;
}

inline
NodeId Graph::num_nodes() const
{
   return _nodes.size();
}

inline
Graph::size_type Graph::num_edges() const
{
   return _num_edges;
}

inline
Node const & Graph::node(NodeId const id) const
{
   // perform index checking
   return _nodes.at(id);
}
//END: Inline section

} // namespace ED

#endif /* GRAPH_HPP */
