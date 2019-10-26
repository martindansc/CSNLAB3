#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <random>

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

// PREPROCESS
Adjlist read_language(string language)
{
    // read the input
    Adjlist adjlist;

    ifstream infile("dependency_networks/" + language);

    string a, b;
    while (infile >> a >> b)
    {
        add_edge(adjlist, a, b);
    }

    return adjlist;
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

// MAIN

void process_language(string language)
{
    Adjlist adjlist = read_language(language);
    Adjlist switching_model = generate_switching_model(adjlist, 5);
    print_network_characteristics(adjlist);
    print_network_characteristics(switching_model);
}

int main()
{
    string language = "Basque_syntactic_dependency_network.txt";

    cout << "Lenguage: " << language << endl;
    process_language(language);
}
