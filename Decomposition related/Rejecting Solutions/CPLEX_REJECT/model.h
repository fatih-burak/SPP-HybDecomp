#ifndef BM_H
#define BM_H

using namespace std;
#include "gurobi_c++.h"
#include "helper_functions.h"
#include "slave_functions.h"
#include <ilcplex/ilocplex.h>

#include <time.h> //Ne pas oublier d'inclure le fichier time.h
#include <iostream> 
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <set>
#include <iterator>     // std::next
#include <map>
#include <list>
#include <tuple>
#include <math.h>
#include <algorithm>


void CPLEX_reject(Instance& inst, Solution& sol);

#endif 

