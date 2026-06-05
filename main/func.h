#ifndef FUNC_H
#define FUNC_H

#include <algorithm>
#include <bitset>
#include <ctime>
#include <iostream>
#include <iterator>
#include <map>
#include <random>
#include <stdexcept>
#include <string>

#include "toylife/toy_plugin.h"
#include "Header.hh"


//std::vector<bool> intToBoolVector(int num, int size);

std::string binaryGenerator(int length = 8, float control = 0.5);

std::map<std::string, int> fromList2Map(std::vector<std::string> list);

std::vector<std::string> decimal2Binary(int numbers);

void printMap(mapa_prot& mapa);
void printMap(mapa_met& mapa);
void printMap(mapa_dim& mapa);
void printMap(mapa_owm& mapa);


std::string fromMapToString(mapa_prot& mapa);
std::string fromMapToString(mapa_met& mapa);
std::string fromMapToString(mapa_dim& mapa);


#endif  // FUNC_H