#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <cmath>
#include <functional>
#include <ctime>

#include "xMath.hpp" // Header file for the function fsolve()
#include "function.hpp"

using namespace std;

int main() {
    clock_t tStart = clock(); // Time
    ofstream altExs("/Users/michaeladmin/Desktop/HRG-AltExS/output/altExshrg.dat");

    vector<vector<double>> data;
    ifstream inputFile("/Users/michaeladmin/Desktop/HRG-AltExS/PDG21Plus_ThFIST.dat");
    if (inputFile.is_open()) {
        string line;
        while (getline(inputFile, line)) {
            vector<double> row;
            stringstream ss(line);
            double value;
            while (ss >> value) {
                row.push_back(value);
            }
            data.push_back(row);
        }
        inputFile.close();
    } else {
        cout << "Failed to open the file." << endl;
        return 1;
    }

    double b = 3.42 * pow(1.0 / 197.3, 3); // MeV
    double a = 329 * pow(1.0 / 197.3, 3); // MeV-2

    for (double T = 10.0; T <= 200; T += 1) {
        for (double muB = 0.0; muB <= 700; muB += 1) {
            double powT3 = pow(T, 3);
                altExs << T << " "<< muB <<" "<< nBB( T,  muB,  b,  a,  data) <<  endl;
               // altExs << T <<" "<< muB << " " << (T/muB)*BardHrg( T,  muB,  b,  a,  data)<< " " << chi2hrgShift( T,  muB,  b,  a,  data)<< endl;

        }
    }

    printf("Time taken: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
    altExs.close();

    return 0;
}
