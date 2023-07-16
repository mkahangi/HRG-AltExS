/**
 * @file main.cpp
 * @brief This code uses Ising model Universality class to introduce a critical point and an alternative expansion scheme to Extend the QCD based Eqaution of state to finite densities and match lattice QCD at low densities
 *
 * @author Micheal Kahangirwe <mkahangi@central.uh.edu>
 *
 * @brief main function
 *
 * Takes inputs from w,  rho  and muBC from the USER
 * The main procedure is:
 *  - Maps Ising co-ordinates  to Alternative expansion scheme then to QCD co-ordinates
 *  - Merges Lattice to Ising
 *  Computes all the Thermodynamic Observables
 *  - Baryon Density, Pressure, Entropy, Baryon number susceptability, Energy, Speed of Sound
 * @return Out puts thermodynamic observables
 **/

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <tuple>

#include "function.hpp"            // Header file for the function fsolve()

//******************************************************************************************************************//
//******************************************************************************************************************//

int main() {
  clock_t tStart = clock(); // Time
  ofstream TrialN1;

  TrialN1.open("/Users/michaeladmin/Desktop/Final_EoS_code/output/TrialN.dat");
  for (double muB = 0; muB <= 700.; muB += 1) {
    for (double T = 5.; T <= 400.; T += 1) {
      TrialN1 << muB << " " << T << " " << Cv(T,muB)<< endl;   
    }
  }
printf("Time taken: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC); 
  TrialN1.close();
  


//******************************************************************************************************************//
// ****************** Checking Stability of the eqaution of state ***************
//******************************************************************************************************************//

  
// Check if the function is stable or unstable
    double minT = 50;    // Minimum value of T
    double maxT = 300;   // Maximum value of T
    double deltaT = 1;  // Step size for T

    double minMuB = 0.0;    // Minimum value of muB
    double maxMuB = 450;    // Maximum value of muB
    double deltaMuB = 1;  // Step size for muB

    bool foundNegative = false;

    // Nested loop to iterate over T and muB values
    for (double T = minT; T <= maxT; T += deltaT) {
        for (double muB = minMuB; muB <= maxMuB; muB += deltaMuB) {
            double result = Cv(T, muB);

            if (result < 0.0) {
                foundNegative = true;
                std::cout << "Negative value found at T = " << T << ", muB = " << muB << " " <<result << std::endl;
                break;  // Exit the inner loop
            }
        }

        if (foundNegative) {
            break;  // Exit the outer loop
        }
    }

    if (!foundNegative) {
        std::cout << "No negative values found for chi2fulln(T, muB)" << std::endl;
    }


    

    return 0;
}







  //******************************************************************************************************************//

