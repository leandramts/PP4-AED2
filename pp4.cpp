#include <iostream>
#include <list>
#include <vector>
#include <stdexcept>
#include <limits>
#include <utility>
#include <string>

using uint = unsigned int;
using Vertex = unsigned int;
using Weight = float;
using VertexWeightPair = std::pair<Vertex, Weight>;

struct Edge {
    Vertex src, dest;
    Weight weight;
};

class WeightedGraphAL
{
private:
    uint num_vertices;
    uint num_edges;
    std::list<VertexWeightPair> *adj;

public:
    WeightedGraphAL(uint num_vertices) : num_vertices(num_vertices),
                                         num_edges(0)
    {
        adj = new std::list<VertexWeightPair>[num_vertices];
    }

    ~WeightedGraphAL()
    {
        delete[] adj;
        adj = nullptr;
    }

    void add_edge(const Vertex &u, const Vertex &v, const Weight &w)
    {
        if (u < 0 || v < 0 || u > num_vertices || v > num_vertices || u == v)
            throw std::invalid_argument("Vertices invalidos");

        auto pair1 = std::make_pair(v, w);

        adj[u].push_back(pair1);

        auto pair2 = std::make_pair(u, w);

        adj[v].push_back(pair2);

        num_edges++;
    }

    const std::list<VertexWeightPair> get_adj(const Vertex &u) const
    {
        if (u < 0 || u > num_vertices)
            throw std::invalid_argument("Vertice invalido");

        return adj[u];
    }

    uint get_num_vertices() const
    {
        return num_vertices;
    }
};

class MinimumPriorityQueue
{
private:
    std::vector<std::pair<Weight, Vertex>> MinHeap;

    int parent(int i)
    {
        return (i - 1) / 2;
    }

    int left(int i)
    {
        return 2 * i + 1;
    }

    int right(int i)
    {
        return 2 * i + 2;
    }

    void min_heapify(int i)
    {
        int l = left(i);
        int r = right(i);
        int smallest = i;

        int tam = MinHeap.size();

        if (l < tam && MinHeap[l].first < MinHeap[i].first) // define qual o menor: o original, esquerda ou direita
            smallest = l;

        if (r < tam && MinHeap[r].first < MinHeap[smallest].first)
        {
            smallest = r;
        }

        if (smallest != i) // se o menor for esquerda ou direita, troca de lugar e recursivamente repete o processo
        {
            std::pair<Weight, Vertex> aux = MinHeap[i];
            MinHeap[i] = MinHeap[smallest];
            MinHeap[smallest] = aux;
            min_heapify(smallest);
        }
    }

public:
    MinimumPriorityQueue() {}

    bool empty()
    {
        return MinHeap.empty();
    }

    std::vector<std::pair<Weight, Vertex>> getHeap() { return MinHeap; }

    void build_min_heap() // constroi a heap binaria minima
    {
        int tam = MinHeap.size();
        for (int i = tam / 2 - 1; i >= 0; i--) // repete o processo de encontrar o menor ate a raiz
        {
            min_heapify(i);
        }
    }

    void insert(std::pair<Weight, Vertex> key)
    {
        MinHeap.push_back(key); // insere chaves na heap
        int i = MinHeap.size() - 1;

        while (i > 0 && MinHeap[i].first < MinHeap[parent(i)].first) // enquanto a chave nao for a raiz e for menor que a pai, eles mudam de lugar
        {
            std::pair<Weight, Vertex> aux = MinHeap[i];
            MinHeap[i] = MinHeap[parent(i)];
            MinHeap[parent(i)] = aux;

            i = parent(i);
        }
    }

    std::pair<Weight, Vertex> extract_min() // remove o elemento de menor chave
    {
        if (MinHeap.size() == 0)
        {
            throw std::runtime_error("Fila vazia");
        }

        std::pair<Weight, Vertex> min = MinHeap[0];
        MinHeap[0] = MinHeap.back();
        MinHeap.pop_back(); // vai pro fim da fila pra ser removido

        if (!MinHeap.empty())
            min_heapify(0);

        return min;
    }
};

class AlgorithmDijkstra
{
private:
    const WeightedGraphAL &g;
    std::vector<uint> dist;
    std::vector<int> pred;

    void initialize(Vertex s)
    {
        uint n = g.get_num_vertices();
        for (Vertex u = 0; u < n; u++)
        {
            if (u == s)
            {
                dist[u] = 0;
                pred[s] = -1;
            }
            else
            {
                dist[u] = std::numeric_limits<uint>::max(); // considera primeiramente distancia infinita
                pred[u] = -1;
            }
        }
    }

    void relax(Vertex u, Vertex v, Weight w)
    {
        if (dist[v] > dist[u] + w)
        {
            dist[v] = dist[u] + w;
            pred[v] = u;
        }
    }

public:
    AlgorithmDijkstra(const WeightedGraphAL &graph) : g(graph)
    {
        dist.resize(g.get_num_vertices());
        pred.resize(g.get_num_vertices());
    }

    void Dijkstra(Vertex s)
    {
        initialize(s);
        MinimumPriorityQueue mq;
        uint n = g.get_num_vertices();

        for (Vertex u = 0; u < n; u++)
        {
            mq.insert({dist[u], u});
        }

        while (!mq.empty())
        {
            auto [d, u] = mq.extract_min();

            for (const auto &adj_vertice : g.get_adj(u))
            {
                Vertex v = adj_vertice.first;
                Weight w = adj_vertice.second;
                if (dist[v] > dist[u] + w)
                {
                    relax(u, v, w);
                    mq.insert({dist[v], v});
                }
            }
        }
    }

    std::vector<Vertex> get_path(Vertex final) const
    {
        std::vector<Vertex> path_inverse;

        if (dist[final] == std::numeric_limits<uint>::max())
        {
            return path_inverse; // quando o vertice final nao e alcancanvel
        }

        Vertex current_vertex = final; // comecamos o processo pelo vertice final

        while (current_vertex != (uint)-1) // a rota do caminho minimo pelo predecessores
        {
            path_inverse.push_back(current_vertex);
            current_vertex = pred[current_vertex];
        }

        std::vector<Vertex> path;

        for (int i = path_inverse.size() - 1; i >= 0; i--) // botar o caminho na ordem correnta
        {
            path.push_back(path_inverse[i]);
        }

        return path;
    }

    const std::vector<uint> &getDistances() const
    {
        return dist;
    }
};

template <typename T>
class UnionFind
{
private:
    std::vector<T> parent;
    std::vector<T> rank;

public:
    UnionFind(T n)
    {
        parent.resize(n);
        rank.resize(n);

        for (T i = 0; i < n; ++i) // make_set() ou initialize()
        {
            parent[i] = i;
            rank[i] = 0;
        }
    }

    void unite(T u, T v) // usando union-by-rank
    {
        T root_u = find(u);
        T root_v = find(v);

        if (root_u != root_v)
        {
            if (rank[root_u] < rank[root_v])
            {
                parent[root_u] = root_v;
            }
            else if (rank[root_u] > rank[root_v])
            {
                parent[root_v] = root_u;
            }
            else
            {
                parent[root_v] = root_u;
                rank[root_u]++;
            }
        }
    }

    T find(T v) // com path compression
    {
        if (v != parent[v])
        {
            parent[v] = find(parent[v]);
        }
        return parent[v];
    }
};

class MSTKruskal {
private:
    unsigned int num_vertices;
    std::vector<Edge> edges;

    void selectionSort(std::vector<Edge>& arr) {
        int n = arr.size();
        for (int i = 0; i < n - 1; i++) {
            int min_idx = i;
            for (int j = i + 1; j < n; j++) {
                if (arr[j].weight < arr[min_idx].weight) {
                    min_idx = j;
                }
            }
            if (min_idx != i) {
                Edge temp = arr[min_idx];
                arr[min_idx] = arr[i];
                arr[i] = temp;
            }
        }
    }

public:
    MSTKruskal(Vertex num_vertices) : num_vertices(num_vertices) {}

    void addEdge(Vertex u, Vertex v, Weight w) {
        edges.push_back({u, v, w});
    }

    std::vector<Edge> findMST() {
        std::vector<Edge> A; 
        UnionFind<Vertex> uf(num_vertices);
        std::vector<Edge> L = edges;
        selectionSort(L);

        for (const auto& edge : L) {
            Vertex root_src = uf.find(edge.src);
            Vertex root_dest = uf.find(edge.dest);
            
            if (root_src != root_dest) {
                A.push_back(edge);
                uf.unite(edge.src, edge.dest);
            }
        }

        return A;
    }
};

class NeuronBlock {
private:
    bool is_diseased;
    MSTKruskal internal_graph;

public:
    NeuronBlock() : is_diseased(false), internal_graph(0) {}

    NeuronBlock(Vertex num_neurons, bool diseased, const std::vector<Edge>& internal_edges) 
            : is_diseased(diseased), internal_graph(num_neurons)
    {
        for (const auto& edge : internal_edges) {
            this->internal_graph.addEdge(edge.src, edge.dest, edge.weight);
        }
    }

    Weight get_sum() {
        std::vector<Edge> mst = internal_graph.findMST();
        Weight sum = 0.0;
        for (const auto& edge : mst) {
            sum += edge.weight;
        }
        return sum;
    }

    bool isDiseased() const {
        return is_diseased;
    }
};

class BrainPaths {
private:
    WeightedGraphAL brain_graph;
    std::vector<NeuronBlock> blocks;
    Vertex num_blocks;

public:
    BrainPaths(Vertex num_brain_blocks) : brain_graph(num_brain_blocks), num_blocks(num_brain_blocks)
    {
        blocks.resize(num_brain_blocks);
    }

    void add_brain_connection(Vertex u, Vertex v, Weight w) {
        brain_graph.add_edge(u, v, w);
    }

    void set_neuron_block(Vertex id, Vertex num_neurons, bool is_diseased, const std::vector<Edge>& edges) {
        blocks[id] = NeuronBlock(num_neurons, is_diseased, edges);
    }

    Weight process_treatment(Vertex entry_point, Vertex exit_point) {
        AlgorithmDijkstra dijkstra(brain_graph);
        dijkstra.Dijkstra(entry_point);
        std::vector<Vertex> path = dijkstra.get_path(exit_point);

        Weight total_weight_sum = 0;

        for (Vertex current_block_id : path) {
                NeuronBlock& current_block = blocks[current_block_id];
                if (current_block.isDiseased()) 
                    total_weight_sum += current_block.get_sum();
        }
            
        return total_weight_sum;
    }
};

int main()
{
    uint num_brain_nodes, num_brain_edges;
    std::cin >> num_brain_nodes >> num_brain_edges;

    BrainPaths brain(num_brain_nodes);

    for (uint i = 0; i < num_brain_edges; i++)
    {
        Vertex u, v;
        Weight w;
        std::cin >> u >> v >> w;
        brain.add_brain_connection(u-1, v-1, w);
    }
    
    Vertex brain_in, brain_out;
    std::cin >> brain_in >> brain_out;

    for (uint i = 0; i < num_brain_nodes; i++)
    {
        uint num_block_nodes, num_block_edges, num_diseased;
        std::cin >> num_block_nodes >> num_block_edges;
        std::cin >> num_diseased;

        if (num_diseased > 0) {
            int dummy_neuron_id;
            for (uint d = 0; d < num_diseased; d++) {
                std::cin >> dummy_neuron_id; // le e descarta os neuronios doentes
            }
        }

        std::vector<Edge> edges;
        for (uint j = 0; j < num_block_edges; j++)
        {
            Vertex u, v;
            Weight w;
            std::cin >> u >> v >> w;
            edges.push_back({u-1, v-1, w});
        }
        brain.set_neuron_block(i, num_block_nodes, num_diseased > 0, edges);
    }

    std::cout << brain.process_treatment(brain_in-1, brain_out-1) << std::endl;
    return 0;
}