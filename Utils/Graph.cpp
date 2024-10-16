#include "Graph.hpp"
#include <limits>
#include <stdexcept>

// Constructor: initializes an adjacency list with a given number of nodes
utils::AdjacencyList::AdjacencyList(size_t nodes_number)
    : nodes(nodes_number)
    , _node_num(nodes_number)
{}

// Function to add a directed edge from one node to another with a specific weight
void utils::AdjacencyList::add_edge(size_t from, size_t to, unsigned weight)
{
    nodes[from].emplace_back(to, weight);
}

// Function to perform topological sorting on the graph
std::vector<size_t> utils::AdjacencyList::topological_sort() const
{
    size_t *indegree = new size_t[_node_num](); // Array to store indegrees of nodes
    for (size_t i = 0; i != _node_num; ++i)
        for (auto edge : nodes[i])
            ++indegree[edge.to]; // Calculate indegree for each node

    std::vector<size_t> topological_sequence; // Vector to store topological order
    topological_sequence.reserve(nodes.size());

    // Iteratively process nodes with zero indegree
    for (size_t pending = _node_num; pending;)
        for (size_t i = 0; i != _node_num; ++i)
            if (indegree[i] == 0)
            {
                indegree[i] = std::numeric_limits<size_t>::max(); // Mark the node as visited
                topological_sequence.emplace_back(i);
                --pending;

                // Decrease indegree of neighboring nodes
                for (auto edge : nodes[i])
                    --indegree[edge.to];
            }

    delete[] indegree; // Free allocated memory
    return topological_sequence;
}

// Function to find the longest path in the graph using topological order
std::vector<size_t> utils::AdjacencyList::get_longest_path() const
{
    std::vector<size_t> path; // Vector to store the longest path
    if (_node_num == 0) return path;

    // Reverse adjacency list to store incoming edges
    std::vector<reverse_node_type> reverse_nodes(_node_num);
    for (size_t from = 0; from != _node_num; ++from)
        for (auto edge : nodes[from])
            reverse_nodes[edge.to].emplace_back(from, edge.weight);

    auto topological_sequence = topological_sort(); // Get topological order of nodes

    auto distance = new size_t[_node_num](); // Array to store distance from source
    auto path_record = new size_t[_node_num]; // Array to record path information

    // Compute longest path using reverse nodes and topological order
    for (size_t to = 1; to != _node_num; ++to)
    {
        const size_t target = topological_sequence[to];
        for (auto edge : reverse_nodes[target])
        {
            unsigned challenger = distance[edge.from] + edge.weight;
            if (challenger > distance[target])
            {
                distance[target] = challenger;
                path_record[target] = edge.from; // Update path
            }
        }
    }

    // Find the end node of the longest path
    size_t end_node = 0;
    for (size_t i = 1; i != _node_num; ++i)
        if (distance[end_node] < distance[i]) end_node = i;

    // Trace back to find the longest path
    for (size_t i = end_node; i != 0; i = path_record[i])
        path.emplace_back(i);

    // Reverse the path to get it in correct order
    for (size_t i = 0, mid = path.size() >> 1; i != mid; ++i)
        std::swap(path[i], path[path.size() - i - 1]);

    delete[] distance;
    delete[] path_record;

    return path;
}

// Function to get the weight of the edge from one node to another
unsigned utils::AdjacencyList::get_weight(size_t from, size_t to) const noexcept
{
    for (auto edge : nodes[from])
        if (edge.to == to) return edge.weight;
}
