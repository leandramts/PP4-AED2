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
        if (u < 0 || v < 0 || u >= num_vertices || v >= num_vertices || u == v)
            throw std::invalid_argument("Vertices invalidos");

        auto pair1 = std::make_pair(v, w);

        adj[u].push_back(pair1);

        auto pair2 = std::make_pair(u, w);

        adj[v].push_back(pair2);

        num_edges++;
    }

    const std::list<VertexWeightPair>& get_adj(const Vertex &u) const
    {
        if (u < 0 || u >= num_vertices)
            throw std::invalid_argument("Vertice invalido");

        return adj[u];
    }

    uint get_num_vertices() const
    {
        return num_vertices;
    }

     uint get_num_edges() const
    {
        return num_edges;
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
    struct Edge {
        Vertex src, dest;
        Weight weight;
    };
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

class Brain
{
    private:
        WeightedGraphAL *g_brain;
        std::vector<WeightedGraphAL*> g_blocks;

        uint order_brain;
        uint tam_brain; 
        std::pair<Vertex, Vertex> in_out_nano_robot;
        std::vector<std::vector<Vertex>> sick_neurons_per_block;
    
    void brain_graph()
    {
        std:: cin >> order_brain >> tam_brain;

        g_brain = new WeightedGraphAL(order_brain);

        for (uint i = 0; i < tam_brain; i++)
        {
            Vertex u, v;
            Weight w;

            std:: cin >> u >> v >> w;
            g_brain->add_edge(u-1, v-1, w);
        }
        
    }

    void init_NanoRobot()
    {
        Vertex in;
        Vertex out;

        std::cin >> in;
        std:: cin >> out;

        in_out_nano_robot = std::make_pair(in-1, out-1);
    }

    void blocks_graph()
    {
        g_blocks.resize(order_brain + 1, nullptr);
        sick_neurons_per_block.resize(order_brain + 1);

        for (uint i = 1; i <= order_brain; i++)
        {
            uint order_block, tam_block;
            std:: cin >> order_block >> tam_block;

            g_blocks[i] = new WeightedGraphAL(order_block);

            uint n_sick_neurons;
            std:: cin >> n_sick_neurons;


            if (n_sick_neurons > 0)
            {
                for (uint s = 0; s < n_sick_neurons; s++)
                {
                    Vertex s_neuron;
                    std::cin >> s_neuron;
                    sick_neurons_per_block[i].push_back(s_neuron - 1);
                }
            }


            for (uint t = 0; t < tam_block; t++)
            {
                Vertex u, v;
                Weight w;
                std:: cin >> u >> v >> w;
                g_blocks[i]->add_edge(u-1, v-1, w);
            }
        }

    }

    public:
        Brain() : g_brain(nullptr) {}
        ~Brain()
    {
        delete g_brain;
        for (WeightedGraphAL* p : g_blocks)
        {
            delete p;
        }
    }


    void input_data()
    {
        brain_graph();
        init_NanoRobot();
        blocks_graph();
    }

};

class NanoRobotsCure
{
    private:
    WeightedGraphAL g;


};

int main()
{
    Brain brain;
    brain.input_data();
    return 0;
}