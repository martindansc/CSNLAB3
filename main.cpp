#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <random>
#include <cmath>
#include <algorithm>
#include <chrono>

using namespace std;

typedef map<string, set<string>> Adjlist;
// MATH

int get_N(Adjlist adjlist)
{
    return adjlist.size();
}

int get_E(Adjlist adjlist)
{
    int sum = 0;
    for (auto const &x : adjlist)
    {
        sum += x.second.size();
    }
    return sum / 2;
}

int get_k(int N, int E)
{
    return 2 * E / N;
}

float get_delta(int N, int E)
{
    return 2 * E / (float)(N * (N - 1));
}

// HELPERS


void print_network_characteristics(Adjlist const &adjlist)
{
    int N = get_N(adjlist);
    int E = get_E(adjlist);
    int k = get_k(N, E);
    float delta = get_delta(N, E);

    cout << "N:" << N << " E: " << E << " <k>: " << k << " delta: " << delta << endl;
}

void print_adjlist(Adjlist adjlist)
{
    for (auto const &x : adjlist)
    {
        cout << x.first << ": ";

        for (auto const &y : x.second)
        {
            cout << y << " ";
        }

        cout << endl;
    }
}

bool add_string_value(Adjlist &adjlist, string word_a, string word_b)
{

    if (adjlist.find(word_a) == adjlist.end())
    {
        adjlist[word_a] = set<string>();
    }

    if (word_a == word_b)
        return false;

    auto const &ret = adjlist[word_a].insert(word_b);

    return ret.second;
}

bool add_edge(Adjlist &adjlist, string a, string b)
{
    add_string_value(adjlist, a, b);
    return add_string_value(adjlist, b, a);
}

bool exists_edge(Adjlist const &adjlist, string a, string b)
{
    auto const &x = adjlist.find(a);
    if (x == adjlist.end())
    {
        return false;
    }

    if (x->second.find(b) == x->second.end())
    {
        return false;
    }

    return true;
}

vector<string> get_node_names(Adjlist adjlist)
{
    vector<string> names;
    for (auto const &x : adjlist)
    {
        names.push_back(x.first);
    }
    return names;
}

bool compare_node_degs(pair<string, int> const &a, pair<string, int> const &b)
{
    return a.second > b.second;
}

bool compare_node_degs_reverse(pair<string, int> const &a, pair<string, int> const &b)
{
    return b.second > a.second;
}

vector<pair<string, int>> get_nodes_with_degs(Adjlist network) 
{
    vector<pair<string, int>> nodes;
    for(auto const &node : network)
    {
        string name = node.first;
        int degree = node.first.size();
        nodes.push_back(pair<string, int>(name, degree));
    }
    return nodes;
}

vector<string> get_strings(vector<pair<string, int>> pairs)
{
    vector<string> tmp;
    for(auto const &pair : pairs)
    {
        tmp.push_back(pair.first);
    }
    return tmp;
}

vector<string> get_node_names_incr(Adjlist network)
{
    auto nodes = get_nodes_with_degs(network);
    sort(nodes.begin(), nodes.end(), &compare_node_degs);
    return get_strings(nodes);
}

vector<string> get_node_names_decr(Adjlist network)
{
    auto nodes = get_nodes_with_degs(network);
    sort(nodes.begin(), nodes.end(), &compare_node_degs_reverse);
    return get_strings(nodes);
}

vector<string> get_node_names_rand(Adjlist network)
{
    auto names = get_node_names(network);
    random_shuffle(names.begin(), names.end());
    return names;
}

// PREPROCESS
Adjlist read_language(string language)
{
    // read the input
    Adjlist adjlist;

    ifstream infile("dependency_networks/" + language);

    string a, b;
    infile >> a >> b; // to ommit the header
    while (infile >> a >> b)
    {
        add_edge(adjlist, a, b);
    }

    return adjlist;
}

// MEASURE
float get_local_clustering(Adjlist adjlist, vector<string> order, bool greater = false, float value = 0.0, float percent = 1.0)
{
    float sum = 0;
    int M = 0;
    int N = get_N(adjlist);
    for (int i = 0; i < order.size()*percent; i++)
    {
        set<string> node = adjlist[order[i]];
        M++;
        int local_sum = 0;
        int num_pairs = 0;
        if (node.size() > 1)
        { // convention that Ci = 0 if degree < 2
            auto y1 = node.begin();
            for (int i = 0; i < node.size(); i++)
            {
                string a = y1->data();

                auto y2 = node.begin();
                advance(y2, i);

                for (int j = i; j < node.size(); j++)
                {
                    string b = y2->data();
                    if (exists_edge(adjlist, a, b))
                    {
                        local_sum++;
                    }
                    num_pairs++;
                    y2++;
                }
                y1++;
            }
            sum += local_sum / (float)num_pairs;
            if (greater)
            {
                if (sum / (float)N + 1 - M / (float)N < value)
                {
                    return 0;
                }
            }
        }
    }
    return sum/(float)N;
}

// SWITCHING MODEL
typedef pair<string, string> Edge;

template <typename Iter, typename RandomGenerator>
Iter select_randomly(Iter start, Iter end, RandomGenerator &g)
{
    std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
    std::advance(start, dis(g));
    return start;
}

template <typename Iter>
Iter select_randomly(Iter start, Iter end)
{
    static std::random_device rd;
    static std::mt19937 gen(rd());
    return select_randomly(start, end, gen);
}

void swap_edge(Edge &edge)
{
    string tmp = edge.first;
    edge.first = edge.second;
    edge.second = tmp;
}

bool is_edge_correct(Edge const &a)
{
    return a.first > a.second;
}

void exchange_ends(Edge &a, Edge &b)
{
    string tmp = a.second;
    a.second = b.second;
    b.second = tmp;
}

void update_edges(Edge &a, Edge &b)
{
    exchange_ends(a, b);
    if (not(is_edge_correct(a)))
        swap_edge(a);
    if (not(is_edge_correct(b)))
        swap_edge(b);
}

vector<Edge>::iterator get_random_edge(vector<Edge> &edges)
{
    return select_randomly(edges.begin(), edges.end());
}

vector<Edge> generate_edge_set(Adjlist const &network)
{
    auto edges = set<Edge>();
    for (auto const &vertex : network)
    {
        string a = vertex.first;

        for (auto const &b : vertex.second)
        {
            auto edge = pair<string, string>(a, b);
            if (not(edge.first > edge.second))
            {
                swap_edge(edge);
            }
            edges.insert(edge);
        }
    }
    return vector<Edge>(edges.begin(), edges.end());
}

bool is_allowed(Adjlist const &network, Edge const &a, Edge const &b)
{
    if (a.first == b.first or a.first == b.second)
        return false;
    else if (a.second == b.first or a.second == b.second)
        return false;
    else if (exists_edge(network, a.first, b.second) or exists_edge(network, b.first, a.second))
        return false;
    else
        return true;
}

void update_network(Adjlist &network, Edge const &a, Edge const &b)
{

    network[a.first].erase(a.second);
    network[a.second].erase(a.first);

    network[b.first].erase(b.second);
    network[b.second].erase(b.first);

    network[a.first].insert(b.second);
    network[b.second].insert(a.first);

    network[b.first].insert(a.second);
    network[a.second].insert(b.first);
}

Adjlist generate_switching_model(Adjlist const &network, unsigned int Q)
{
    Adjlist switching_model = Adjlist(network);
    vector<Edge> edges = generate_edge_set(switching_model);
    unsigned int E = edges.size();
    for (int i = 0; i < Q * E; i++)
    {
        Edge &a = *get_random_edge(edges);
        Edge &b = *get_random_edge(edges);

        if (is_allowed(switching_model, a, b))
        {
            update_network(switching_model, a, b);
            update_edges(a, b);
        }
    }
    return switching_model;
}
// GRAPH GENERATION

Adjlist generate_erdos_renyi(vector<string> names, int N, int E)
{
    Adjlist adjlist;
    vector<string> aux;

    int i = 0;
    while (i < E)
    {
        int a = rand() % N;
        int b = rand() % N;

        if (add_edge(adjlist, names[a], names[b]))
        {
            i++;
        }
    }
    return adjlist;
}

// ALGORITHMS
float p_monte_carlo_erdos_renyi(Adjlist adjlist, float x, int N, int E, int T)
{
    vector<string> names = get_node_names(adjlist);
    int times = 0;
    for (int i = 0; i < T; i++)
    {
        Adjlist random_graph = generate_erdos_renyi(names, N, E);
        float c_value = get_local_clustering(random_graph, names);
        if (c_value > x)
            times++;
    }

    return times / (float)T;
}

float p_monte_carlo_switching(Adjlist adjlist, float x, int N, int E, int T)
{
    vector<string> names = get_node_names(adjlist);
    int times = 0;
    for (int i = 0; i < T; i++)
    {
        Adjlist random_graph = generate_switching_model(adjlist, log(E));
        float c_value = get_local_clustering(random_graph, names);
        if (c_value > x)
            times++;
    }

    return times / (float)T;
}

float p_monte_carlo_switching_exact_optimization(Adjlist adjlist, float x, int N, int E, int T)
{
    vector<string> names = get_node_names(adjlist);
    int times = 0;
    for (int i = 0; i < T; i++)
    {
        Adjlist random_graph = generate_switching_model(adjlist, log(E));
        float c_value = get_local_clustering(random_graph, names, true, x);
        if (c_value > x)
            times++;
    }

    return times / (float)T;
}

float p_monte_carlo_switching_approximate_optimization(Adjlist adjlist, float x, int N, int E, int T)
{
    vector<string> names = get_node_names(adjlist);
    int times = 0;
    for (int i = 0; i < T; i++)
    {
        Adjlist random_graph = generate_switching_model(adjlist, log(E));
        float c_value = get_local_clustering(random_graph, names, true, x, 0.1);
        if (c_value > x)
            times++;
    }

    return times / (float)T;
}


float p_monte_carlo_erdos_renyi_exact_optimization(Adjlist adjlist, float x, int N, int E, int T)
{
    vector<string> names = get_node_names_decr(adjlist);
    int times = 0;
    for (int i = 0; i < T; i++)
    {
        Adjlist random_graph = generate_erdos_renyi(names, N, E);
        float c_value = get_local_clustering(random_graph, names, true, x);
        if (c_value > x)
            times++;
    }

    return times / (float)T;
}


// MAIN

void process_language(string language)
{
    Adjlist adjlist = read_language(language);
    Adjlist switching_model = generate_switching_model(adjlist, 20);
    float x = get_local_clustering(adjlist, get_node_names(adjlist));
    int N = get_N(adjlist);
    int E = get_E(adjlist);
    int k = get_k(N, E);
    float delta = get_delta(N, E);
    
    chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    p_monte_carlo_switching_exact_optimization(adjlist, x, N, E, 1);
    chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    cout << "Time difference = " << chrono::duration_cast<chrono::milliseconds> (end - begin).count() << "[ms]" << endl;
}

void perform_ordering_tests(string const &language, vector<string> (*ordering)(Adjlist))
{
    Adjlist adjlist = read_language(language);
    float x = get_local_clustering(adjlist, get_node_names_decr(adjlist));
    int N = get_N(adjlist);
    int E = get_E(adjlist);
    int k = get_k(N, E);
    float delta = get_delta(N, E);
    Adjlist switching = generate_switching_model(adjlist, 20);
    Adjlist erdos = generate_erdos_renyi(get_node_names(adjlist), N, E);
    vector<string> order = ordering(switching);
    {
        cout << "Erdos: " << endl;
        chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        get_local_clustering(erdos, order, true, x);
        chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        cout << "Time difference = " << chrono::duration_cast<chrono::milliseconds> (end - begin).count() << "[ms]" << endl;
    }
    {
        cout << "Switching: " << endl;
        chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        get_local_clustering(switching, order, true, x);
        chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        cout << "Time difference = " << chrono::duration_cast<chrono::milliseconds> (end - begin).count() << "[ms]" << endl;
    }
}

int main()
{
    std::vector<std::string> languages = {"English.txt", "Greek.txt", "Hungarian.txt", "Italian.txt", "Turkish.txt"};
    std::vector<vector<string>(*)(Adjlist)> orderings;
    orderings.push_back(get_node_names);
    orderings.push_back(get_node_names_rand);
    orderings.push_back(get_node_names_incr);
    orderings.push_back(get_node_names_decr);
    for (auto const &language : languages)
    {
        cout << "language: " << language << endl;
        int i = 0;
        for (auto const &ordering : orderings)
        {
            cout << "Ordering " << i++ << endl;
            perform_ordering_tests(language, ordering);
        }
        
    }

}
