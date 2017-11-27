#ifndef READ_GRAPH_HPP
#define READ_GRAPH_HPP

#include "graph.hpp"

#include <utility> //std::pair

namespace ED {
std::pair<int, ED::Graph> read_graph(std::string filename);
}

#endif
