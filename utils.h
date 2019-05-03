#ifndef utils_H
#define utils_H

#include "particle.h"
#include "vose.h"

#include <vector>
#include <string>
#include <cmath>
#include <cstdio>
#include <chrono>

#include <iostream>

//Define a function to check if a directory
//exists for the write state function
#include <sys/types.h>
#include <sys/stat.h>

using namespace std;

bool dir_check(char const *pathname);

vector<particle> initialize_system(int part_num, double mass, double length);

void update_state(vector<particle> &particles, const int Ng, const double dt, const double length,
                  double *charge, double *potential, double *EField);
void write_state(vector<particle>&, double*, double *, double*, int, int i);

double temperature(vector<particle>&);
void scale_mom(vector<particle>&, double);

void calc_charge(int Ng, double dx, double *charge, vector<particle> particles);
void calc_potential(int Ng, double dx, double *potential, double *charge);
void calc_EField(int Ng, double dx, double *EField, double *potential);

#endif //utils_H
