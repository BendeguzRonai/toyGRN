/*
* Based on Prokaryotes by Samuel von der Dunk
* https://github.com/samvonderdunk/Prokaryotes
*/
#include "Prokaryote.hh"
#include "stdio.h"
#include <stdlib.h>
#include <png.h>
#include <zlib.h>



Prokaryote::Prokaryote() {
}

Prokaryote::~Prokaryote() {
	if (G != NULL) {
        delete G; 
        G = NULL;
    }
}

void Prokaryote::EmptyProkaryote()
{
    G = NULL;
	fossil_id = 0;
	time_of_appearance = 0;
	Ancestor = NULL;
	mutant = false;
	alive = true;
	agent = NULL;
	//fitness assigned later
}

void Prokaryote::InitialiseProkaryote(){
	EmptyProkaryote();
	mutant = true;	//The first will be deemed a mutant so that the ancestor trace always comes back to at least one individual of the initial batch.
	if(genome_init != "")	 G = new string(genome_init);
	else	{
		string rand_g = "";
		for (int i = 0; i < genome_length; i++) {
			//float numero = (float)rand() / (float)RAND_MAX;
			
			rand_g +=  (uniform() > 0.5) ? '1' : '0';
			
		}
        G = new std::string(rand_g);
	}
}

void Prokaryote::Mutate() {
	
    for (size_t i = 0; i < G->length(); ++i) {

        if (uniform() < point_mut) {
            (*G)[i] = ((*G)[i] == '1') ? '0' : '1'; 
			mutant = true;
			relative_distance = 1.;
			fitness = 0.;
        }
    }
	if (!chrom) return;
	if (uniform() < del_mut && G->length()>40){
		size_t start_idx = static_cast<size_t>(uniform() * (G->length() - 19)); // Random valid start index
		(*G).erase(start_idx, 20); // Remove 20 characters
		mutant = true;
		relative_distance = 1.;
		fitness = 0.;
	}
	if (uniform() < copy_mut && G->length()<200){
		size_t start_idx = static_cast<size_t>(uniform() * (G->length() - 19));
        std::string copied_segment = G->substr(start_idx, 20);
        *G += copied_segment; // Append the copied segment
        mutant = true;
		relative_distance = 1.;
        fitness = 0.;
	}
	if (uniform() < cut_mut){
		size_t start_idx = static_cast<size_t>(uniform() * (G->length() - 19));
		std::string cut_segment = G->substr(start_idx, 20);
		*G += cut_segment;
		(*G).erase(start_idx, 20);
		mutant = true;
		relative_distance = 1.;
		fitness = 0.;
	}
}

void Prokaryote::FitCalc(){

	agent = std::make_shared<toyAgent>(G);
	agent->RegulatoryPhenotype(toy);
	int e = agent->LogicAssert(); //currently only for 1st 2 genes
	if (e>0) { 
		std::cout << "Logic Function Error: " << e << endl;
		agent->Print();
	}
	
	if (Ancestor!=NULL && agent->logic_function == Ancestor->agent->logic_function && agent->morphogen_outleak==Ancestor->agent->morphogen_outleak)
	{	
		agent = Ancestor->agent;
		relative_distance = Ancestor->relative_distance;
		fitness = Ancestor->fitness;
	}
	else {
		agent->PatternFormation();
		relative_distance = agent->CalculateRelativeDistance();
		fitness = agent->FitnessCalculation();
		//agent->Print();
	}
	//PrintData(true);
	// delete &agent;
	// agent=NULL;
	// printf("Logic function:\n");

	// for(it = LogicFunction.begin(); it != LogicFunction.end(); it++, i++) {
	// 	// Increment distance for each differing character
	// 	printf("%d ", LogicFunction[i]);
	// }
	// printf("##############################################\n");

	
	// int H_dist = 0;
	// vector<int>::iterator it;
	// int i = 0;

	// for(it = LogicFunction.begin(); it != LogicFunction.end(); it++, i++) {
	// 	// Increment distance for each differing character
	// 	if (LogicFunction[i] != target[i]) H_dist++;
	// }

	// double f = 1.0 - static_cast<double>(H_dist) / (*G).length();
	// fitness = f*f;

}

// void Prokaryote::FitCalc(){
	
// 	int H_dist = 0;

// 	for (size_t i = 0; i < (*G).length(); ++i) {
//         // Increment distance for each differing character
//         if ((*G)[i] != target[i]) {
//             H_dist++;
//         }
//     }

// 	double f = 1.0 - static_cast<double>(H_dist) / (*G).length();
// 	fitness = f*f;
// }


void Prokaryote::Clone(Prokaryote* parent, unsigned long long tot_prok_count)
{
	EmptyProkaryote();

	G = new string(*parent->G);
	
	fossil_id = tot_prok_count;
	time_of_appearance = Time;
	relative_distance = parent->relative_distance;	
	fitness = parent->fitness;
	agent = parent->agent;
	
	if (parent->mutant)	Ancestor = parent;	//If your parent was a mutant (i.e. its genome holds a mutation with respect to its parent), then your immediate ancestor is your parent.
	else	Ancestor = parent->Ancestor;	//If your parent was not a mutant (its genome is the same as its parent), then point to its MRCA.
}

void Prokaryote::PrintData(bool include_agent_data)
{
	printf("##############################################\n");
	printf("Prokaryote #%llu:\n", fossil_id);
	if(Ancestor==NULL)	printf("Generation 0\n");
	else	printf("Child of #%llu\n", Ancestor->fossil_id);
	printf("Time of birth = %d\n", time_of_appearance);
	printf("Relative distance = %f\n", relative_distance);
	printf("Fitness = %f\n", fitness);
	printf("It is %s.\n", (alive)? "alive":"dead");
	printf("It is %sa mutant.\n", ((mutant)? "":"not " ));
	printf("Genome:\n%s\n", G->c_str());
	
	
	if (agent!=NULL && include_agent_data) agent->Print();

	
	printf("##############################################\n");
}
