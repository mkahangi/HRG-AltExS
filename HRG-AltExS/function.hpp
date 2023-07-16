#ifndef FUNCTION_HPP
#define FUNCTION_HPP

#pragma once

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <tuple>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>



using namespace std;


double PHIn(double T, double mi, double di);

double PHIBaryon(double T, const vector<vector<double>>& data);

double nBB(double T, double muB, double b, double a, const vector<vector<double>>& data);

double nAB(double T, double muB, double b, double a, const vector<vector<double>>& data);

double BardHrg(double T, double muB, double b, double a, const vector<vector<double>>& data);

double chi2HrgA(double T, double muB, double b, double a, const vector<vector<double>>& data);

double dchi2HrgdTA(double T, double muB, double b, double a, const vector<vector<double>>& data);

double d2chi2HrgdT2A(double T, double muB, double b, double a, const vector<vector<double>>& data);

double chi3HrgA(double T, double muB, double b, double a, const vector<vector<double>>& data);

double chi4HrgA(double T, double muB, double b, double a, const vector<vector<double>>& data);

double chi5HrgA(double T, double muB, double b, double a, const vector<vector<double>>& data);


double chi6HrgA(double T, double muB, double b, double a, const vector<vector<double>>& data);

double Kappa2(double T, double muB, double b, double a, const vector<vector<double>>& data);

double Kappa4(double T, double muB, double b, double a, const vector<vector<double>>& data);

double chi2hrgShift(double T, double muB, double b, double a, const vector<vector<double>>& data)


#endif
