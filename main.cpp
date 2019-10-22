#include <fstream>
#include <string>
#include <iostream>

using namespace std;

void process_language(string language) {
    ifstream infile("dependency_networks/" + language);

    string a, b;
    while (infile >> a >> b)
    {
        cout << a << " " << b << endl;
    }
}

int main() {
    string language = "Basque_syntactic_dependency_network.txt";
    process_language(language);
}


