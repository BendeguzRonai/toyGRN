/*
* Based on Prokaryotes by Samuel von der Dunk
* https://github.com/samvonderdunk/Prokaryotes
*/
#ifndef GeneralHeader
#define GeneralHeader

using namespace std;
//using namespace __gnu_cxx;

#include <stdio.h>
#include <stdlib.h>
#include <ext/numeric>
#include <cmath>
#include <list>
#include <vector>
#include <set>
#include <iostream>
#include <iomanip>
#include <map>
#include <iterator>
#include <algorithm>
#include <boost/next_prior.hpp>
#include <cstring>
#include <sstream>
#include <fstream>
#include "dSFMT.h"	// Kirsten used it for some functions, so lets just stick with it
#include <sys/time.h>
#include <ctime>
#include <typeinfo>
#include <string>

#include "toylife/helper_functions.h"
#include "toylife/toy_plugin.h"
// #include "/home/sam/Programmes/nvwa-1.1/nvwa/debug_new.h"

#define toDigit(c) (c-'0')  // Converts char to digit

//SIMULATION PARAMETERS
const int NR=50;  
const int NC=50;
const int TimeZero=0;
const int default_SimTime=1000;
const int TimeTerminalOutput = 1000;  
const int TimeSaveGrid = 1000; 
const int TimePruneFossils = 1000;
const int TimeOutputFossils = 1000000;
//const bool quickExp = false;

//DNA & Mutations 
const double point_mut=0.0001;
const double del_mut=0.0005;
const double copy_mut=0.0005;
const double cut_mut=0.0005;
extern bool chrom; // chormosomal mutations

const bool DiverseInitialPop=false; //Each cell gets its own random genome. false for uniform population w/ 1 random genome
const int genome_length = 40; //Must be multiple of 20 
extern string genome_init;

//Patterns
const int tmax = 100; //time
const int ncells=31; //ncells
extern vector<vector<vector<int>>> Patterns;
//const double base_fitness = 0.15;

//constants for Population.cc
const double diffusion_rate = 0.;  // >1: multiple Margolus steps per time step, <1: probability of single Margolus step each time step.
const double death_rate = 0.1;
const int replication_neighbourhood = 3; 


extern double selection_strength; // parameter used in the fitness calculation function

//Variables defined in main.cc
extern int Time;
extern int initial_seed;
extern unsigned long long seed_draws;
extern string folder;
extern string backup_reboot;
extern string anctrace_reboot;
extern int SimTime;
extern int Target;
extern bool well_mixing;

//Mode: Genotype space sampling
//extern bool sim;
extern std::string mode;
extern int N;

extern dsfmt_t dsfmt;
inline double uniform()
{
  seed_draws ++;
  return dsfmt_genrand_close_open(&dsfmt);
}

const string genome_file="";
const string genestate_file="";
const string backup_file="";
const string anctrace_file="";
extern const ToyPlugin toy;


const bool vocal=false;



#endif
