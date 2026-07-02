/*
* Based on Prokaryotes by Samuel von der Dunk
* https://github.com/samvonderdunk/Prokaryotes
*/
#ifndef ProkaryoteHeader
#define ProkaryoteHeader

#include "Header.hh"
#include "stdlib.h"
#include "toylife/helper_functions.h"
#include "toylife/toy_plugin.h"
#include "func.h"
#include "toyAgent.hh"



class Prokaryote{
	public:
		string* G;
		std::shared_ptr<toyAgent> agent;
		
		double relative_distance;
		double fitness;

		//For the fossil record.
		int time_of_appearance;
		unsigned long long fossil_id;	//Now it should be 64-bit (32-bit unsigned would already be about 200 times as big as the last id you get out of a 100x100 simulation of 100k AUT).
		Prokaryote* Ancestor; //
		bool mutant;
		bool alive;
		typedef std::list<Prokaryote*>::iterator iterpps;

		Prokaryote();
		~Prokaryote();

		void EmptyProkaryote();
		void InitialiseProkaryote();
		void Clone(Prokaryote* PP_template, unsigned long long tot_prok_count);
		void Mutate();
		void FitCalc();
		void PrintData(bool include_genome_data);

		

};

#endif
