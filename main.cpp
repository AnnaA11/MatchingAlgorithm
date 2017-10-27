#include <iostream>
#include <cstdlib>

#include "graph.hpp"

//! return graph consisting of path with \c num_nodes many vertices
static ED::Graph create_path(ED::NodeId num_nodes)
{
   ED::Graph result{num_nodes};

   for (ED::NodeId node_id = 0; node_id + 1 < num_nodes; ++node_id)
   {
      result.add_edge(node_id, node_id + 1);
   }
   // will be moved
   return result;
}

int main(int argc, char**)
{
   // create and output a path with argc many vertices
   ED::Graph const graph = create_path(static_cast<ED::NodeId>(argc));
   std::cout << graph;
   return EXIT_SUCCESS;
}
