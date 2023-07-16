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
#include "function.hpp" // Header file for the function fsolve()

using namespace std;

double PHIn(double T, double mi, double di) {
    const double PI = 3.14159265359;
    return (di * pow(mi, 2) * T) / (2 * PI * PI) * xMath::BesselK(2, mi / T);
}

double PHIBaryon(double T, const vector<vector<double>>& data) {
    double result = 0.0;
    int dataSize = data.size();

    for (int i = 0; i < dataSize; i++) {
         //cout <<" "<<data[i][1] <<" "<< " "<<data[i][2]<< " "<<data[i][4] << endl;
        if (data[i][4] == 1) {
            result += PHIn(T, data[i][1] * 1000, data[i][2]);
           
        }
    }
    return result;
}



//double nBB(double T, double muB, double b, double a, const vector<vector<double>>& data) {
//    double nBr = 50000.0;
//    double tolerance = 1e-8;  // Tolerance for convergence
//    int maxIterations = 100000;  // Maximum number of iterations
//
//    // Define the function to find the root of
//    auto func = [&](double nBr) {
//        return (b * nBr * exp((b * nBr) / (1 - b * nBr) - (2 * a * nBr) / T)) / (1 - nBr * b) - b * PHIBaryon(T, data) * exp(muB / T);
//    };
//
//    // Apply the Secant algorithm
//    int iter = 0;
//    double nBr_prev = nBr + 1.0;  // Initial guess for the previous iteration
//
//    while (iter < maxIterations) {
//        double f = func(nBr);
//        double dfdx = (f - func(nBr_prev)) / (nBr - nBr_prev);
//        double delta = -f / dfdx;
//        nBr_prev = nBr;
//        nBr += delta;
//
//        // Check convergence
//        if (abs(delta) < tolerance)
//            break;
//
//        ++iter;
//    }
//
//    return nBr;
//}

double nBB(double T, double muB) {
    const double a = 1.0; // Define the value of 'a'
    const double b = 1.0; // Define the value of 'b'
    const double epsilon = 1e-6; // Define the desired accuracy of the solution

    // Initial guess for nBr
    double nBr = 50000.0;

    // Newton-Raphson iteration
    while (true) {
        double numerator = b * nBr * exp((b * nBr) / (1 - b * nBr) - (2 * a * nBr) / T);
        double denominator = 1 - nBr * b;
        double lhs = numerator / denominator;
        double rhs = b * PHIBaryon(T) * exp(muB / T);

        double f = lhs - rhs;
        double df = (numerator / denominator + (b * (b * nBr) * exp((b * nBr) / (1 - b * nBr) - (2 * a * nBr) / T)) / pow(1 - b * nBr, 2)) / T;

        double delta = f / df;
        nBr -= delta;

        // Check convergence
        if (std::abs(delta) < epsilon) {
            break;
        }
    }

    return nBr;
}

double nAB(double T, double muB, double b, double a, const vector<vector<double>>& data) {
    double nB = 50000;

    double tolerance = 1e-8;  // FindRoot tolerance
    int maxIterations = 1000;  // maximum number of iterations

    // Define the function to find the root of
    auto func = [&](double nB) {
        return (b * nB * exp((b * nB) / (1 - b * nB) - (2 * a * nB) / T)) / (1 - nB * b) - b * PHIBaryon(T, data) * exp(-muB / T);
    };

    // Apply the FindRoot algorithm
    for (int iter = 0; iter < maxIterations; ++iter) {
        double f = func(nB);
        double dfdx = (func(nB + tolerance) - func(nB - tolerance)) / tolerance;
        double delta = f / dfdx;
        nB -= delta;

        // Check convergence
        if (abs(delta) < tolerance)
            break;
    }
    return nB;
}
double BardHrg(double T, double muB, double b, double a, const vector<vector<double>>& data){
    return nBB( T,  muB,  b,  a, data) - nAB( T,  muB,  b,  a,  data);
}
/// Analytic derivatives of baryon density////////
double chi2HrgA(double T, double muB, double b, double a, const vector<vector<double>>& data){
     double nBBN = nBB( T,  muB,  b,  a,  data);
    return (T*nBBN*pow(-1 + b*nBBN,2))/(T - 2*a*nBBN*pow(-1 + b*nBBN,2));
}
double dchi2HrgdTA(double T, double muB, double b, double a, const vector<vector<double>>& data) {
     double diff=0.001;
    return (chi2HrgA(T- 2 * diff, muB, b, a, data) - 8 * chi2HrgA(T- diff, muB, b, a, data) +
                8 * chi2HrgA(T + diff, muB,b, a, data) - chi2HrgA(T + 2 * diff, muB, b, a, data)) / (12 * diff);
}

double d2chi2HrgdT2A(double T, double muB, double b, double a, const vector<vector<double>>& data) {
     double diff=0.001;
    return (2*chi2HrgA(T- 3 * diff, muB, b, a, data) - 27*chi2HrgA(T- 2 * diff, muB, b, a, data) +
    270*chi2HrgA(T- 1 * diff, muB, b, a, data) - 490*chi2HrgA(T, muB, b, a, data) +
    270*chi2HrgA(T+ 1* diff, muB, b, a, data) - 27*chi2HrgA(T+ 2 * diff, muB, b, a, data) +
    2*chi2HrgA(T+ 3 * diff, muB, b, a, data))/(180*pow(diff,2));
    }

double chi3HrgA(double T, double muB, double b, double a, const vector<vector<double>>& data){
     double nBBN = nBB( T,  0,  b,  a,  data);
    return -((pow(T,3)*nBBN*(1 - 3*b*nBBN)*pow(-1 + b*nBBN,3))/pow(T - 2*a*nBBN*pow(-1 + b*nBBN,2),3));
}

double chi4HrgA(double T, double muB, double b, double a, const vector<vector<double>>& data){
     double nBBN = nBB( T,  0,  b,  a,  data);
    return (pow(T,4)*nBBN*pow(-1 + b*nBBN,4)*(3*T*pow(1 - 3*b*nBBN,2) -
       2*(1 - 4*b*nBBN + 6*pow(b,2)*pow(nBBN,2))*(T - 2*a*nBBN*pow(-1 + b*nBBN,2))))/
   pow(T - 2*a*nBBN*pow(-1 + b*nBBN,2),5);
}

double chi5HrgA(double T, double muB, double b, double a, const vector<vector<double>>& data){
     double nBBN = nBB( T,  0,  b,  a,  data);
    return (T*pow(-1 + b*nBBN,3)*(-12*pow(T,5)*nBBN*pow(-1 + b*nBBN,2)*(-1 + 3*b*nBBN)*(1 - 4*b*nBBN + 6*pow(b,2)*pow(nBBN,2)) +
       6*pow(T,4)*nBBN*pow(-1 + b*nBBN,2)*(T - 2*a*nBBN + 4*a*b*pow(nBBN,2) - 2*a*pow(b,2)*pow(nBBN,3))*
        (-1 + 5*b*nBBN - 10*pow(b,2)*pow(nBBN,2) + 10*pow(b,3)*pow(nBBN,3)) -
       ((1 - 3*b*nBBN)*(3*pow(T,6)*nBBN*pow(1 - 3*b*nBBN,2)*pow(-1 + b*nBBN,2) +
            4*pow(T,4)*(T*nBBN*pow(-1 + b*nBBN,2))/(T - 2*a*nBBN*pow(-1 + b*nBBN,2))*(T - 2*a*nBBN*pow(-1 + b*nBBN,2))*
             (3*T*pow(1 - 3*b*nBBN,2) - 2*(1 - 4*b*nBBN + 6*pow(b,2)*pow(nBBN,2))*(T - 2*a*nBBN*pow(-1 + b*nBBN,2)))))/
        (T - 2*a*nBBN*pow(-1 + b*nBBN,2))))/pow(T - 2*a*nBBN*pow(-1 + b*nBBN,2),6);
}
double chi6HrgA(double T, double muB, double b, double a, const vector<vector<double>>& data){
     double nBBN = nBB( T,  0,  b,  a,  data);
     double nBB2N = chi2HrgA( T,  0,  b,  a,  data);
     double nBB3N = chi3HrgA( T,  0,  b,  a,  data);
     double nBB4N = chi4HrgA( T,  0,  b,  a,  data);
     double nBB5N = chi5HrgA( T,  0,  b,  a,  data);

    return -(nBB2N*((10*nBB3N*nBB4N*(1 - 3*b*nBBN))/(pow(nBBN,2)*pow(-1 + b*nBBN,3)) +
       (5*nBB2N*nBB5N*(1 - 3*b*nBBN))/(pow(nBBN,2)*pow(-1 + b*nBBN,3)) +
       (10*nBB2N*(3*pow(nBB3N,2) + 2*nBB2N*nBB4N)*(1 - 4*b*nBBN + 6*pow(b,2)*pow(nBBN,2)))/
        (pow(nBBN,3)*pow(-1 + b*nBBN,4)) + 24*pow(nBB2N,5)*
        (pow(nBBN,-5) + (5*pow(b,6)*nBBN)/pow(-1 + b*nBBN,6) - (6*pow(b,5))/pow(-1 + b*nBBN,5)) +
       10*pow(nBB2N,3)*nBB3N*(-6/pow(nBBN,4) - (24*pow(b,5)*nBBN)/pow(-1 + b*nBBN,5) +
          (30*pow(b,4))/pow(-1 + b*nBBN,4))));
          }
double Kappa2(double T, double muB, double b, double a, const vector<vector<double>>& data){
     return chi4HrgA(T, muB, b, a, data)/(6*T*dchi2HrgdTA(T, muB, b, a,data));
}
double Kappa4(double T, double muB, double b, double a, const vector<vector<double>>& data){
     return 1/(360*T*pow(dchi2HrgdTA(T, muB, b, a,data),3))*(3*pow(dchi2HrgdTA(T, muB, b, a, data),2)*chi6HrgA(T, muB, b, a, data) - 5*d2chi2HrgdT2A(T, muB, b, a, data)*pow(chi4HrgA(T, muB, b, a, data),2));
}
double chi2hrgShift(double T, double muB, double b, double a, const vector<vector<double>>& data){
    return 2*chi2HrgA(T*(1 + Kappa2( T,  0,  b,  a, data)*pow(muB/T,2) + Kappa4( T,  0,  b,  a, data)*pow(muB/T,4) ), 0, b, a, data);
}

