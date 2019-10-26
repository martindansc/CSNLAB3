#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <set>
#include <map>

using namespace std;

typedef map<string, set<string> > Adjlist;

// HELPERS

void print_adjlist(Adjlist adjlist) {
    for (auto const& x : adjlist) {
        cout << x.first << ": ";

        for(auto const& y : x.second) {
            cout << y << " ";
        }

        cout << endl;
    }
}

bool add_string_value(Adjlist& adjlist, string word_a, string word_b) {

    if(adjlist.find(word_a) == adjlist.end()) {
        adjlist[word_a] = set<string>();
    }

    if(word_a == word_b) return false;

    auto const &ret = adjlist[word_a].insert(word_b);

    return ret.second;
}

bool add_edge(Adjlist& adjlist, string a, string b) {
    add_string_value(adjlist, a, b);
    return add_string_value(adjlist, b, a);   
}

bool exists_edge(Adjlist& adjlist, string a, string b) {
    auto const &x = adjlist.find(a);
    if(x == adjlist.end()) {
        return false;
    }
    
   if(x->second.find(b) == x->second.end()) {
       return false;
   }
   
   return true;
}

vector<string> get_node_names(Adjlist adjlist) {
    vector<string> names;
    for (auto const& x : adjlist) {
        names.push_back(x.first);
    }
    return names;
}

// MATH

int get_N(Adjlist adjlist) {
    return adjlist.size();
}

int get_E(Adjlist adjlist) {
    int sum = 0;
    for (auto const& x : adjlist) {
        sum += x.second.size();
    }
    return sum/2;
}

int get_k(int N, int E) {
    return 2*E/N;
}

float get_delta(int N, int E) {
    return 2*E/(float)(N*(N-1));
}


// PREPROCESS
Adjlist read_language(string language) {
    // read the input
    Adjlist adjlist;

    ifstream infile("dependency_networks/" + language);

    string a, b;
    while (infile >> a >> b) {
        add_edge(adjlist, a, b);
    }

    return adjlist;
}


// GRAPH GENERATION

Adjlist generate_erdos_renyi(vector<string> names, int N, int E) {
    Adjlist adjlist;

    vector<string> aux;

    int i = 0;
    while(i < E) {
        int a = rand() % N;
        int b = rand() % N;

        if(add_edge(adjlist, names[a], names[b])) {
            i++;
        }
    }

    return adjlist;
}

// MAIN

void process_language(string language) {
    Adjlist adjlist = read_language(language);

    int N = get_N(adjlist);
    int E = get_E(adjlist);
    int k = get_k(N, E);
    float delta = get_delta(N, E);

    cout << "N:" << N << " E: " << E << " <k>: " << k << " delta: " << delta << endl; 

    // compute......
    cout << exists_edge(adjlist, ",", "azaldu") << endl;
    cout << exists_edge(adjlist, ",", "xxxxx") << endl;

    Adjlist random_graph = generate_erdos_renyi(get_node_names(adjlist), N, E);

    print_adjlist(random_graph);
}

int main() {
    string language = "Basque_syntactic_dependency_network.txt";

    cout << "Lenguage: " << language << endl;
    process_language(language);
}

