/*
* Based on Prokaryotes by Samuel von der Dunk
* https://github.com/samvonderdunk/Prokaryotes
*/
#include "FossilRecord.hh"

FossilRecord::FossilRecord()
{
}

FossilRecord::~FossilRecord()
{
	iterpps ips;
	ips = FossilList.begin();
	
	while (ips != FossilList.end())
	{	
		if (*ips != NULL)
		{	
			delete (*ips);
			ips = FossilList.erase(ips);
		}
		else	ips++;
	}
}

void FossilRecord::EraseFossil(unsigned long long fossilID)
{
	iterpps ips;
	ips = FossilList.begin();
	while (ips != FossilList.end())
	{
		if ((*ips)->fossil_id == fossilID)
		{
			ips = FossilList.erase(ips);
			return;
		}
		ips++;
	}
}

void FossilRecord::BuryFossil(Prokaryote *P)
{
	FossilList.push_back(P);
}

void FossilRecord::ExhibitFossils()
{
	FILE* f;
	char OutputFile[800];
	sprintf(OutputFile, "%s/ancestors/anctrace%08d.txt", folder.c_str(), Time);
	f=fopen(OutputFile, "w");

	if (f == NULL)	printf("Failed to open file for writing the ancestor trace.\n");
	fprintf(f, "#id\t#anc_id\t#time_oa\t#fit\t#genome\n");	

	iterpps ip = FossilList.begin();
	while(ip != FossilList.end())
	{
		if ((*ip)->Ancestor == NULL)
		{
			fprintf(f, "%llu\t%d\t%d\t%f\t%s\n", (*ip)->fossil_id, 0, (*ip)->time_of_appearance, (*ip)->fitness, (*ip)->G->c_str());
			// fprintf(f, "----------Expression Pattern----------\n");
			// for (size_t t = 0; t < (*ip)->agent->Pattern[0].size(); ++t) {  // Loop over columns (100 time steps)
			// 	for (size_t n = 0; n < (*ip)->agent->Pattern.size(); ++n) fprintf(f, "%d", (*ip)->agent->Pattern[n][t]);  
			// 	fprintf(f, "\n");  
			// }
			fprintf(f, "----------Phenotypic Pattern----------\n");

			for (size_t t = 0; t < (*ip)->agent->PhenotypicPattern[0].size(); ++t) {  // Loop over columns (100 time steps)
				for (size_t n = 0; n < (*ip)->agent->PhenotypicPattern.size(); ++n) fprintf(f, "%d", (*ip)->agent->PhenotypicPattern[n][t]);  
				fprintf(f, "\n");  
			}
		}
		else
		{
			fprintf(f, "%llu\t%llu\t%d\t%f\t%s\n", (*ip)->fossil_id, ((*ip)->Ancestor)->fossil_id, (*ip)->time_of_appearance, (*ip)->fitness, (*ip)->G->c_str());
			
			// fprintf(f, "----------Expression Pattern----------\n");
			// for (size_t t = 0; t < (*ip)->agent->Pattern[0].size(); ++t) {  // Loop over columns (100 time steps)
			// 	for (size_t n = 0; n < (*ip)->agent->Pattern.size(); ++n) fprintf(f, "%d", (*ip)->agent->Pattern[n][t]);  
			// 	fprintf(f, "\n");  
			// }
			fprintf(f, "----------Phenotypic Pattern----------\n");

			for (size_t t = 0; t < (*ip)->agent->PhenotypicPattern[0].size(); ++t) {  // Loop over columns (100 time steps)
				for (size_t n = 0; n < (*ip)->agent->PhenotypicPattern.size(); ++n) fprintf(f, "%d", (*ip)->agent->PhenotypicPattern[n][t]);  
				fprintf(f, "\n");  
			}
		}
		ip++;
	}
	fclose(f);
}



				

