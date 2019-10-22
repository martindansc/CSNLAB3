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

int get_delta(int N, int E) {
    return 2*E/(N*(N-1));
}


// PREPROCESS

void add_string_value(Adjlist& adjlist, string word_a, string word_b) {

    if(adjlist.find(word_a) == adjlist.end()) {
        cout << word_a << endl;
        adjlist[word_a] = set<string>();
    }

    if(word_a == word_b) return;

    adjlist[word_a].insert(word_b);
}

Adjlist read_language(string language) {
    // read the input
    Adjlist adjlist;

    ifstream infile("dependency_networks/" + language);

    string a, b;
    while (infile >> a >> b) {
        add_string_value(adjlist, a, b);
        add_string_value(adjlist, b, a);   
    }

    return adjlist;
}


// MAIN

void process_language(string language) {
    Adjlist adjlist = read_language(language);

    // compute......
}

int main() {
    string language = "Basque_syntactic_dependency_network.txt";
    process_language(language);
}

