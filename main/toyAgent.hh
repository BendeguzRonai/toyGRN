#ifndef Agent
#define Agent

#include <algorithm>  // std::count std::min
#include <bitset>     // std::bitset
#include <iostream>
#include <stdexcept>  // std::invalid_argument
#include <string>
#include <vector>
#include "Header.hh"
#include "stdlib.h"

// toyLife code
#include "toylife/helper_functions.h"
#include "toylife/toy_plugin.h"
#include "func.h"

class toyAgent {
   public:

    std::string* genotype;

    // toyLife things
    mapa_prot prots;
    mapa_dim dims;
    mapa_met mets;
    mapa_owm owns;

    vector<vector<int>> logic_function;
    vector<int> morphogen_outleak;


    std::vector<std::pair<int, int> > vec_genotype;
    int g;
    int state_space;
    std::vector<std::vector<int>> Pattern;
    std::vector<std::vector<int>> PhenotypicPattern;
    

    toyAgent(string* genotype);
    ~toyAgent();

    void Print();

    void reacting(const ToyPlugin& toy);

    int promoter_expression(const ToyPlugin& toy);

    vector<std::string> getAgentData();

    //COMPUTE PHENOTYPE, moved from toy plugin

    void RegulatoryPhenotype(const ToyPlugin& toy);
    double FitnessCalculation();
    int CalculateHammingDistance();
    double CalculateRelativeDistance();

    void PatternFormation();
    void PrintPattern();

    int LogicAssert();
    int FunctionToInt();
    
    // bool finding_cycles(
    //     vector<pair<int, pair<mapa_prot, mapa_dim> > >&
    //         iterations,
    //     int expression);

    // vector<pair<int, pair<mapa_prot, mapa_dim> > > boolean_network(const ToyPlugin& toy);



};

#endif