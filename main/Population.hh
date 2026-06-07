/*
* Based on Prokaryotes by Samuel von der Dunk
* https://github.com/samvonderdunk/Prokaryotes
*/
#ifndef PopulationHeader
#define PopulationHeader

#include "Prokaryote.hh"
#include "Header.hh"
#include "FossilRecord.hh"
#include <cstdio>
#include <stdlib.h>

class Population
{
	public:
		Prokaryote* PPSpace[NR][NC];
		FossilRecord* Fossils;

		int p_nr_proks_;	//Probably need this...

		typedef std::list<int>::iterator iter;
		typedef std::list<unsigned long long>::iterator iterull;
		typedef std::list<Prokaryote*>::iterator iterpps;
		typedef std::pair<int,int> coords;	//Allow to define row, col pair of int's (use for functions).

		unsigned long long p_id_count_;	// Counter for all agents

		//These are for looking at the occurrence of evolution.
		
		double MeanDistance;
		double LeastDistance;
		double MeanFitness;	//Add noise for a given period of time.
		double MaxFitness;
		int size;
		Population();
		~Population();

        void InitialisePopulation();
        void Reset();
		
		void ContinuePopulationFromBackup();
		void ReadBackupFile();
		void ReadAncestorFile();

		void Stats();


		void UpdatePopulation();
		void MargolusDiffusion();
		void WellMix();


		coords PickNeighbour(int i, int j);
		int NeighbourhoodDensity(int i, int j);


		void PruneFossilRecord();

		void PrintFieldToFile();
		
};

#endif
