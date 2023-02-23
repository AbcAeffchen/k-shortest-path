#ifndef KSP_GRAPHRW_H
#define KSP_GRAPHRW_H

#endif //KSP_GRAPHRW_H

#include "Graph.h"
#include "GraphGenerator.h"

#include <cassert>
#include <cstring>
#include <fstream>
#include <sstream>

/**
 * Read and write graph files.
 */
class GraphRW
{
public:
//    enum class WeightConversion {KEEP_ORIGINAL, REDRAW_UNIFORM, REDRAW_EXPONENTIAL, REDRAW_POWERLAW};
    enum class Weights : int {keep, Uniform, Random, SumDeg, ProdDeg, InvSumDeg, InvProdDeg};

private:
    template<bool prepareForDeltaStepping>
    static constexpr NodeType getID(NodeType node)
    {
        if constexpr(prepareForDeltaStepping)
            return 2 * node;
        else
            return node;
    }

    template<EdgeConcept EdgeType>
    static void preparingForDeltaStepping(std::vector<EdgeIdType>& outEdgeBeginnings, std::vector<EdgeType>& edges, const typename EdgeType::WeightType delta)
    {
        const NodeType numNodes = (static_cast<NodeType>(outEdgeBeginnings.size()) - 1) / 2;
        for(NodeType i = 0; i < numNodes; i++)
        {
            const auto firstHeavyEdge = std::partition(edges.begin() + static_cast<EdgeIdDifferenceType>(outEdgeBeginnings[getID<true>(i)]),
                                                       edges.begin() + static_cast<EdgeIdDifferenceType>(outEdgeBeginnings[getID<true>(i + 1)]),
                                                       [delta](const auto& edge) { return edge.getWeight() < delta; });

            std::sort(edges.begin() + static_cast<EdgeIdDifferenceType>(outEdgeBeginnings[getID<true>(i)]),
                      edges.begin() + static_cast<EdgeIdDifferenceType>(outEdgeBeginnings[getID<true>(i + 1)]));

            outEdgeBeginnings[getID<true>(i) + 1] = static_cast<EdgeIdType>(std::distance(edges.begin(), firstHeavyEdge));
        }
    }

public:

    static void writeMetisFile(const bool directed, const bool weighted,
                               const NodeType numNodes, const EdgeIdType numEdges,
                               const GraphGenerator::AdjacencyListType& adjacencyList,
                               const std::string& path, const unsigned long seed = 0)
    {
        // write to file
        std::ofstream metisFile;
        metisFile.open(path);

        metisFile << "% " << (directed ? "directed" : "undirected") << "\n";

        if(seed > 0)
            metisFile << "% Seed: " << seed << "\n";

        metisFile << numNodes << " " << (directed ? numEdges : numEdges / 2) << " " << (weighted ? 1 : 0) << "\n"; // the 1 is the flag for weighted graphs

        for(const auto& neighbors : adjacencyList)
        {
            for(const auto neighbor : neighbors)
            {
                // Metis stores nodes indexed at 1
                metisFile << (neighbor.getTarget() + 1) << " ";
                if(weighted)
                    metisFile << neighbor.getWeight() << " ";
            }

            metisFile << "\n";
        }

//        metisFile << std::endl;

        metisFile.close();
    }

private:
    [[nodiscard]] static bool fileTypeIs(const std::string& file, const std::string& type) noexcept
    {
        if(file.length() >= type.length())
            return 0 == file.compare(file.length() - type.length(), type.length(), type);
        else
            return false;
    }

public:
    template<typename WT, bool directed, bool prepareForDeltaStepping = false, Weights WC = Weights::keep>
    static auto readFile(const std::string& file, WT delta = static_cast<WT>(0.1))
    {
        if(fileTypeIs(file, ".metis"))
            return readMetisFile<WT, directed, prepareForDeltaStepping, WC>(file, delta);
        else if(fileTypeIs(file, ".gr"))
            return readDIMACSFile<WT, directed, prepareForDeltaStepping, WC>(file, delta);

        Print::info() << file << " has no supported format." << std::endl;
        exit(1);
    }

    template<typename WT, bool directed, bool prepareForDeltaStepping = false, Weights WC = Weights::keep>
    static auto readDIMACSFile(const std::string& file, WT delta = static_cast<WT>(0.1))
    {
        /*
         * A graph contains n nodes and m arcs
         * Nodes are identified by integers 1...n
         * Graphs can be interpreted as directed or undirected, depending on the problem being studied
         * Graphs can have parallel arcs and self-loops
         * Arc weights are signed integers
         * http://www.diag.uniroma1.it//challenge9/format.shtml
         *  c This is a comment
         *  p sp n m
         *  a U V W
         */

        std::ifstream graphFile(file, std::ifstream::in);

        graphFile.seekg(0);     // move to the beginning of the file
        std::string line;

        // get problem line
        do
        {
            std::getline(graphFile, line);
        }
        while(line.at(0) == 'c');

        std::istringstream fileStream(line);
        std::string tmp, lineType;
        NodeType numNodes;
        EdgeIdType numEdges;
        // read first line with the settings "p sp numNodes numEdges"
        fileStream >> lineType >> tmp >> numNodes >> numEdges;

        using EdgeType = Edge<double, true>;
        using AdjacencyListType = std::vector<std::vector<EdgeType>>;

        AdjacencyListType neighbors(numNodes);

        size_t minWeight = std::numeric_limits<size_t>::max();
        size_t maxWeight = 0;

        while(std::getline(graphFile, line))
        {
            // skip empty lines and comment lines starting with a "c"
            if(line.length() == 0 || line.at(0) == 'c')
                continue;

            std::stringstream ss(line);
            NodeType u, v;
            size_t weight;  // DIMACS stores only integer weights
            minWeight = std::min(minWeight, weight);
            maxWeight = std::max(maxWeight, weight);

            ss >> lineType;
            if(lineType != "a")
            {
                Print::info() << file << " has invalid format." << std::endl;
                exit(1);
            }

            ss >> u >> v >> weight;
            // node indexes need to be adjusted to start at 0.
            neighbors[u - 1].emplace_back(v - 1, static_cast<double>(weight));
        }

        // convert weights to floats in [0,1]
        for(auto& uNeighbors : neighbors)
        {
            for(auto& e : uNeighbors)
            {
                e.setWeight(e.getWeight() / static_cast<double>(maxWeight));
            }
        }

        std::vector<size_t> outEdgeBeginnings(prepareForDeltaStepping ? (2 * numNodes + 1) : (numNodes + 1), 0);
        std::vector<EdgeType> edges;
        edges.reserve(numEdges);

        for(NodeType u = 0; u < neighbors.size(); u++)
        {
            outEdgeBeginnings[2 * u] = edges.size();
            for(const auto& e : neighbors[u])
            {
                edges.push_back(e);
            }
        }
        outEdgeBeginnings.back() = edges.size();

        if constexpr(prepareForDeltaStepping)
        {
            preparingForDeltaStepping(outEdgeBeginnings, edges, delta);
        }

        if constexpr(prepareForDeltaStepping)
            return Graph<Edge<WT, directed>, true>(numNodes, outEdgeBeginnings, edges, static_cast<WT>(delta), 1.0);
        else
            return Graph<Edge<WT, directed>, false>(numNodes, outEdgeBeginnings, edges);
    }

    static auto convertEdgelistFile(const std::string& file, const Weights weightType = Weights::keep, const bool makeUndirected = true)
    {
        /*
         * comment line starts with '%'
         * other lines contain source-target pairs
         * nodes numbered beginning at 1.
         */

        std::ifstream graphFile(file, std::ifstream::in);

        graphFile.seekg(0);     // move to the beginning of the file
        std::string line;

        using EdgeType = Edge<float, true>;
        using AdjacencyListType = std::vector<std::vector<EdgeType>>;

        std::vector<std::tuple<NodeType, NodeType, float>> edgelist;

        NodeType maxNodeId = 0;

        while(std::getline(graphFile, line))
        {
            // skip empty lines and comment lines starting with a "c"
            if(line.length() == 0 || line.at(0) == '%')
                continue;

            std::stringstream ss(line);
            NodeType u, v;
            ss >> u >> v;

            maxNodeId = std::max(maxNodeId, v - 1);

            edgelist.template emplace_back(u - 1, v - 1, 1.0);
            if(makeUndirected)
                edgelist.template emplace_back(v - 1, u - 1, 1.0);
        }

        std::ranges::sort(edgelist);

        maxNodeId = std::max(maxNodeId, std::get<0>(edgelist.back()));
        const NodeType numNodes = maxNodeId + 1;
        AdjacencyListType neighbors(numNodes);

        // convert weights to floats in [0,1]
        for(const auto& [u, v, weight] : edgelist)
        {
            neighbors[u].emplace_back(v, weight);
        }

        std::random_device seedGen;
        std::default_random_engine gen(seedGen());
        std::uniform_real_distribution<float> drawWeight(0, 1);

        // set weights
        float maxWeight = 0;
        for(NodeType u = 0; u < neighbors.size(); u++)
        {
            const auto uDeg = static_cast<NodeType>(neighbors[u].size());

            for(size_t i = 0; i < neighbors[u].size(); i++)
            {
                const auto vDeg = static_cast<NodeType>(neighbors[neighbors[u][i].getTarget()].size());
                float weight = 1;

                switch(weightType)
                {
                    case Weights::Uniform:
                        // nothing to do
                        break;
                    case Weights::keep:
                        // nothing to do
                        break;
                    case Weights::Random:
                        weight = drawWeight(gen);
                        break;
                    case Weights::SumDeg:
                        weight = static_cast<float>(uDeg + vDeg);
                        break;
                    case Weights::ProdDeg:
                        weight = static_cast<float>(uDeg * vDeg);
                        break;
                    case Weights::InvSumDeg:
                        weight = 1.0f / static_cast<float>(uDeg + vDeg);
                        break;
                    case Weights::InvProdDeg:
                        weight = uDeg * vDeg == 0 ? 1.0f : 1.0f / static_cast<float>(uDeg * vDeg);
                        break;
                }

                maxWeight = std::max(maxWeight, weight);

                neighbors[u][i].setWeight(weight);
            }
        }

        // normalize weights
        if(maxWeight > 1)
        {
            for(NodeType u = 0; u < neighbors.size(); u++)
            {
                for(size_t i = 0; i < neighbors[u].size(); i++)
                {
                    neighbors[u][i].setWeight(neighbors[u][i].getWeight() / maxWeight);
                }
            }
        }

        GraphGenerator::getGiantComponentDirected(neighbors);

        return neighbors;
    }

    template<typename WT, bool directed, bool prepareForDeltaStepping = false, Weights WC = Weights::keep>
    static auto readMetisFile(const std::string& file, WT delta = static_cast<WT>(0.1))
    {
        static_assert(std::is_same<WT, void>::value || WeightConcept<WT>);
        // WT == void ? int : WT
        constexpr bool weighted = !std::is_same<WT, void>::value;
        using WeightType = typename std::conditional<weighted, WT, int>::type;

        /*
         * Note about METIS:
         * The first line contains the number of nodes, the number of edges and flag
         * The flag can be 0 (= unweighted) or 1 (= weighted)
         * The i-th (not commented) line contain a list of the neighbors of the
         * i-th node. Notice that the nodes are numbered from 1 to n (not 0 to n-1).
         * The number of edges does not include backward edges, so this number has
         * to be doubled.
         */

        std::ifstream graphFile(file, std::ifstream::in);

        graphFile.seekg(0);     // move to the beginning of the file
        std::string line;

        // get the head of the file containing the number of nodes and edges and a flag
        // skip comment lines at the beginning
        do
        {
            std::getline(graphFile, line);
        }
        while(line.at(0) == '%');

        std::istringstream fileStream(line);
        NodeType numNodes;
        EdgeIdType numEdges;
        int fileType;
        // read first line with the settings
        fileStream >> numNodes >> numEdges >> fileType;

        if(fileType > 1)
            throw std::invalid_argument("Unknown METIS file type flag.");

        if(weighted && fileType == 0)
            std::cerr << "The file is unweighted but is stored as weighted graph. This needs more memory as needed." << std::endl;

        if(!directed)
            numEdges *= 2;      // undirected edges are stored as two directed edges

        // initialize lists
        std::vector<EdgeIdType> outEdgeBeginnings;
        outEdgeBeginnings.resize(getID<prepareForDeltaStepping>(numNodes) + 1);         // adding a sentinel

        std::vector<Edge<WT, directed>> edges;
        edges.reserve(numEdges);

        NodeType src = 0;
        WeightType heaviestWeight = 0;

        while(std::getline(graphFile, line) && src < numNodes)
        {
            // An empty line represents a node with no outgoing edges. Continue with the next node.
            if(line.length() == 0)
            {
                outEdgeBeginnings[getID<prepareForDeltaStepping>(src)] = edges.size();
                src++;
                continue;
            }

            // skip comment lines
            if(line.at(0) == '%')
                continue;

            // set the current node
            outEdgeBeginnings[getID<prepareForDeltaStepping>(src)] = edges.size();

            NodeType dest;
            auto weight = static_cast<WeightType>(1);      // if weighted it gets overwritten, if not, every weight is either set to 1 or not stored at all
            std::stringstream ss(line);

            while(ss >> dest)
            {
                if(fileType == 1)
                {
                    ss >> weight;
                }

                if constexpr(weighted)
                {
                    if(WC != Weights::keep)
                    {
                        assert(false);      // todo not implemented yet
                    }

                    if constexpr(prepareForDeltaStepping)
                        heaviestWeight = std::max(heaviestWeight, weight);

                    edges.emplace_back(dest - 1, weight);     // METIS labels nodes from 1 to n, we count nodes from 0 to n - 1
                }
                else
                {
                    edges.emplace_back(dest - 1);
                }
            }

            src++;
        }

        assert(getID<prepareForDeltaStepping>(numNodes) == outEdgeBeginnings.size() - 1);
        assert(numEdges == edges.size() || numEdges * 2 == edges.size());

        // setting the sentinel
        outEdgeBeginnings[getID<prepareForDeltaStepping>(numNodes)] = edges.size();

        if(directed && numEdges * 2 == edges.size())
        {
            Print::info() << std::endl << "WARNING: Undirected graph read as directed." << std::endl;
        }

        if constexpr(prepareForDeltaStepping)
        {
            preparingForDeltaStepping(outEdgeBeginnings, edges, delta);
        }

        // link all backward edges to the corresponding forward edges.
        if(!directed)
        {
            assert(false);
            // todo implement this
            // GraphRW::link_backward_edges(nodes, edges);
        }

        if constexpr(prepareForDeltaStepping)
            return Graph<Edge<WT, directed>, true>(numNodes, outEdgeBeginnings, edges,
                                                   static_cast<WT>(delta), heaviestWeight);
        else
            return Graph<Edge<WT, directed>, false>(numNodes, outEdgeBeginnings, edges);
    }
};