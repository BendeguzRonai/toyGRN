/*
* Based on Prokaryotes by Samuel von der Dunk
* https://github.com/samvonderdunk/Prokaryotes
*/
#ifndef FossilHeader
#define FossilHeader

#include "Prokaryote.hh"
#include "Header.hh"
#include <cstdio>
#include <stdlib.h>

class FossilRecord
{
	public:
		std::list<Prokaryote*> FossilList;

		typedef std::list<Prokaryote*>::iterator iterpps;
  	FossilRecord();
  	~FossilRecord();
		void EraseFossil(unsigned long long fossilID);
  	void BuryFossil(Prokaryote *P);
		void ExhibitFossils();
};
#endif
