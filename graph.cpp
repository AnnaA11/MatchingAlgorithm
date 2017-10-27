#include "graph.hpp" // always include corresponding header first

#include <ostream>
#include <stdexcept>

namespace ED
{
/////////////////////////////////////////////
//! \c Node definitions
/////////////////////////////////////////////

void Node::add_neighbor(NodeId const id)
{
   _neighbors.push_back(id);
}

/////////////////////////////////////////////
//! \c Graph definitions
/////////////////////////////////////////////

Graph::Graph(NodeId const num_nodes)
   :
   _nodes(num_nodes),
   _num_edges(0)
{}

void Graph::add_edge(NodeId node1_id, NodeId node2_id)
{
   if (node1_id == node2_id)
   {
      throw std::runtime_error("ED::Graph class does not support loops!");
   }

   // minimum redundancy :-), maybe a bit overkill...
   auto impl = [this](NodeId a, NodeId b)
   {
      Node & node = _nodes.at(a);
      node.add_neighbor(b);
   };

   impl(node1_id, node2_id);
   impl(node2_id, node1_id);

   ++_num_edges;
}

std::ostream & operator<<(std::ostream & str, Graph const & graph)
{
   str << "c This encodes a graph in DIMACS format\n"
       << "p edge " << graph.num_nodes() << " " << graph.num_edges() << "\n";

   for (NodeId node_id = 0; node_id < graph.num_nodes(); ++node_id)
   {
      auto const & node = graph.node(node_id);

      for (auto const & neighbor_id : node.neighbors())
      {
         // output each edge only once
         if (node_id < neighbor_id)
         {
            str << "e " << to_dimacs_id(node_id) << " " << to_dimacs_id(neighbor_id) << "\n";
         }
      }
   }

   str << std::flush;
   return str;
}


/////////////////////////////////////////////
//! global functions
/////////////////////////////////////////////

NodeId from_dimacs_id(DimacsId const dimacs_id)
{
   if (dimacs_id == 0)
   {
      throw std::runtime_error("Invalid (0) DIMACS id.");
   }

   return static_cast<NodeId>(dimacs_id - 1);
}

DimacsId to_dimacs_id(NodeId const node_id)
{
   if (node_id == std::numeric_limits<NodeId>::max())
   {
      throw std::runtime_error("Invalid (inf) node id.");
   }

   return static_cast<DimacsId>(node_id + 1);
}

} // namespace ED
