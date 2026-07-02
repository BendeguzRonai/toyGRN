/*
* Based on Prokaryotes by Samuel von der Dunk
* https://github.com/samvonderdunk/Prokaryotes
*/
#include "Header.hh"
#include "Population.hh"
#include "Prokaryote.hh"
#include "dSFMT.h"
#include "toylife/toy_plugin.h"
#include "toyAgent.hh"


dsfmt_t dsfmt;
int Time;
int initial_seed = time(0);
bool chrom = false;
bool PrintGenome = false;
string gnm = "";
std::string genome_init = "";
string genomes[1000000];

unsigned long long seed_draws = 0;
string folder = "../../Results/";
string OutputFile = "out_dir/GenomeData";

std::string mode = "standard";
int N=1000000;
int Target = 0;
bool well_mixing = false;
int SimTime = default_SimTime;

double selection_strength = 1; // parameter used in the fitness calculation function


string genestate_init = genestate_file;
string backup_reboot = backup_file;
string anctrace_reboot = anctrace_file;

void Setup(int argc, char** argv);

const ToyPlugin toy;  

vector<vector<vector<int>>> Patterns;
vector<int> Rules(165536, -1); //Rules[i] == PatternID??
vector<int> PatternAbundance(42638, 0); //in logic space. PatterAbundance[paternID]== #rules


void readPatterns(string allpatterns, string rulefile) {

	ifstream infile(allpatterns.c_str());
	ifstream infile2(rulefile.c_str());
	string line;
	char*  data;
	int ipatt, npatt, success;
	int i = 0;
	vector<vector<int>> Pattern;

	if (!infile.is_open())
	{
		printf("Pattern file could not be opened.\n");
		exit(1);
	}


	while(getline(infile,line))
	{
		if (i%101==0){

			data = (char*)line.c_str();
			success = sscanf(data, "%d %d\n", &ipatt, &npatt);

			if(success != 2 or i/101!=ipatt)
			{
				cerr << "Could not read pattern ID\n" << endl;
				exit(1);
			}

			if (i!=0) Patterns.push_back(Pattern);
			Pattern.clear();

			PatternAbundance[ipatt] = npatt;
		}
		else{
			stringstream ss(line);
            vector<int> row;
            int val;

            while (ss >> val) {
                row.push_back(val);
            }
            Pattern.push_back(row);
        }
		i++;
	}


	if (!infile2.is_open())
	{
		printf("Rule file could not be opened.\n");
		exit(1);
	}
    while (getline(infile2, line)) {

        stringstream ss(line);
        vector<int> ruleSet;
        int num;

        while (ss >> num) {
            ruleSet.push_back(num);
        }
        
		for (size_t i = 1; i < ruleSet.size(); i++) {
        Rules[ruleSet[i]] = ruleSet[0];
		}
    }
    //for (size_t i = 0; i < Rules.size(); i++) cout << Rules[i] << " ";
    infile.close();
	infile2.close();
}

void printParsedData(const vector<vector<vector<int>>>& patterns, const vector<int>& ruleCounts) {
    cout << "\nParsed file1.txt (Patterns and Rule Counts):\n";
    for (size_t i = 0; i < 10; i++) {
        cout << "Pattern " << i << " (Rules: " << ruleCounts[i] << ")\n";
        for (const auto& row : patterns[i]) {
            for (int val : row) {
                cout << val << " ";
            }
            cout << "\n";
        }
        cout << "\n";
    }
}

void SampleGenomes(int N, string filename){
	/*
	1. Read the patterns
	2. Generate N random genomes length 40
	3. Find the logic rule & make pattern
	4. Check if pattern matches file
	5. If so, keep count of patterns discovered + #new per timestep 
	*/
	Prokaryote* PP;
	int PPfunction;

	vector<int> DiscoveredPatternAbundance(42638, 0);
	vector<int> DiscoveryRate(N/1000+1, 0);
	int PhenoCounter = 0;
	
	for (int i=1; i<N+1; i++){
		if (i%10000 == 0) cout << i << " "<< PhenoCounter<< "\n";
		if (i%1000==0) DiscoveryRate[i/1000] = PhenoCounter;
		PP = new Prokaryote();
		PP->InitialiseProkaryote();
		PP->agent = std::make_shared<toyAgent>(PP->G);
		PP->agent->RegulatoryPhenotype(toy);
		PP->agent->PatternFormation();

		PPfunction = PP->agent->FunctionToInt();

		//Compare patterns
		int H=0;
		for (int t = 0; t < tmax; ++t) {
			for (int k = 0; k < ncells; ++k) {
				if (Patterns[Rules[PPfunction]][t][k] != PP->agent->Pattern[k][t]) {
					++H;
					cout << k << "," << t <<" " << Patterns[Rules[PPfunction]][t][k] << " " << PP->agent->Pattern[k][t] << " ";
				}  
			}
		}

		if (H>0){
			printf("Assertion error\n");
			cout << PPfunction << " " << Rules[PPfunction] << "\n";
			printf("PP\n");
			
			PP->agent->Print();

			printf("Pattern in dataset\n");
			for (const auto& row : Patterns[Rules[PPfunction]]) {
				for (int val : row) {
					cout << val << " ";
				}
				cout << "\n";
			}
			cout << "Dist: "<< H << "\n";
			//break;
		}
		else{
			++DiscoveredPatternAbundance[Rules[PPfunction]];
			if (DiscoveredPatternAbundance[Rules[PPfunction]]==1) ++PhenoCounter;
			
		}
		delete PP;
		PP = NULL;
	} 

	std::ofstream file(filename);
	if (!file) printf("No fiiileeeee\n");  // Check if file opened successfully

	for (size_t i = 0; i < DiscoveredPatternAbundance.size(); ++i) {
		if (DiscoveredPatternAbundance[i]>0) file << i << "\t" << DiscoveredPatternAbundance[i] << "\n";  // Write each number followed by a space
	}
	
	for (int num : DiscoveryRate) {
		file << num << " ";  // Write each number followed by a space
	}
}

int main(int argc, char** argv) {
	
	printf("\n\033[93m### Setup ###\033[0m\n");
	Setup(argc, argv);
	dsfmt_init_gen_rand(&dsfmt, initial_seed);	//Used to seed uniform().
	srand(initial_seed);	//Used to seed random_shuffle(...).

	printf("\b\nSetup completed...\n\n");

	if (Target!=0){
		printf("\n\033[93m### Reading Patterns ###\033[0m\n");
		readPatterns("PatternsData/PatternsLeak0.txt", "PatternsData/RulesLeak0.txt");
		
	}

    //######################################################################################
    //############################## SAMPLE GENOMES MODE ###################################
    //######################################################################################
	if (mode=="GPsample"){
		
		string filename = folder + "/GPscan" + to_string(N/1000000) +"M" +to_string(initial_seed) +".txt";
		SampleGenomes(N, filename);
		return 0;
	}
    //######################################################################################
    //############################## PRINT GENOMES MODE ####################################
    //######################################################################################
	else if (mode=="printGenome"){
		

		//Just 1 genome
		if (gnm != "" && gnm.find('.') == std::string::npos) {
			
			std::shared_ptr<toyAgent> agent = std::make_shared<toyAgent>(&gnm);
			agent->RegulatoryPhenotype(toy);
			agent->PatternFormation();
			if (Target!=0){
                double f = agent->FitnessCalculation();
			    agent->Print();
			    cout << "Fitness / target: " << f << " / " << Target << "\n"; 
            } 
            else cout << "You forgot to give me a Target pattern. If you dont I crash :(\n"; 
            
			return 0;
		}

		//Multiple genomes. Have to decide what to do with this
		else if (gnm.find('.') != std::string::npos){
			std::ifstream file(gnm);
            if (!file) {
                std::cerr << "Failed to open genome file: " << gnm << std::endl;
            }
            std::string line;
			int i =0;
            while (std::getline(file, line)) {
                genomes[i] = line; // std::getline removes the newline character
				i++;
            }
            file.close();
		}
		
		
		FILE* f=fopen(OutputFile.c_str(), "w");
		fprintf(f, "Gnm\tGenes\tOutlieak\tLogic\n");

		for (int i=0; i<1000000; i++)
		{	
			vector<int>::iterator it;
			int j;

			if (genomes[i] != "")
			{	
				//cout << genomes[i] << "\n";
				std::shared_ptr<toyAgent> agent = std::make_shared<toyAgent>(&genomes[i]);
				agent->RegulatoryPhenotype(toy);
				
				if (Target!=0){
					agent->PatternFormation();
					double fit = agent->FitnessCalculation();
					cout << fit << ","; 
				}

				else {

					fprintf(f, "%s\t", genomes[i].c_str());
					for (int p = 0; p < agent->vec_genotype.size(); ++p) fprintf(f, "%d,",toy.prot_gen[agent->vec_genotype[p].second]);
					fprintf(f, "\t");
					j=0;
					for(it = agent->morphogen_outleak.begin(); it != agent->morphogen_outleak.end(); it++, j++) fprintf(f, "%d,", agent->morphogen_outleak[j]);
					fprintf(f, "\t");

					for (int leak_in=0; leak_in < 3; leak_in++){
						j=0;
						for(it = agent->logic_function[leak_in].begin(); it != agent->logic_function[leak_in].end(); it++, j++) fprintf(f, "%d,", agent->logic_function[leak_in][j]);
					}
					fprintf(f, "\n");
				}


			}
		}
		cout << "\n";
		return 0;
	}
    //######################################################################################
    //########################### MULTIPLE EXPERIMENTS MODE ################################
    //######################################################################################
    else if (mode=="quickSim"){
        Population* P;	


        for(int experiment=0; experiment < 20; experiment++){

            cout << "Experiment " << experiment << "\n";
            
        
            P = new Population();
            P->InitialisePopulation();

            for(Time=0; Time<SimTime+1; Time++){	
                if(Time%TimePruneFossils==0 && Time!=0)	P->PruneFossilRecord();
                if (Time%100==0){
                    P->Stats();
                    
					if(P->MaxFitness==1) {
						cout << "Time: " << Time <<" Mean fitness: " << P->MeanFitness << " Max fitness: " << P->MaxFitness << " Mean distance: " << P->MeanDistance << " Least distance: " << P->LeastDistance << " Size: " << P->size << endl;
						break;
                    }
                }
                if (Time%10000==0){
                    P->Stats();
                    cout << "Time: " << Time <<" Mean fitness: " << P->MeanFitness << " Max fitness: " << P->MaxFitness << " Mean distance: " << P->MeanDistance << " Least distance: " << P->LeastDistance << " Size: " << P->size << endl;
                 } 
            
                if (P->size==0) {
                    cout << "Time: " << Time <<" Mean fitness: " << P->MeanFitness << " Max fitness: " << P->MaxFitness << " Mean distance: " << P->MeanDistance << " Least distance: " << P->LeastDistance << " Size: " << P->size << endl;
                    break;
                }
                P->UpdatePopulation();		//Maybe move to main?????
                
            }    
            P->PruneFossilRecord();
            //P->Stats();
			//cout << "Time: " << Time <<" Mean fitness: " << P->MeanFitness << " Max fitness: " << P->MaxFitness << " Size: " << P->size << endl;
            delete P;
	        P = NULL;
        }
        return 0;
    }
    //######################################################################################
    //########################### COMPARE PHENOTYPES FROM 2 GENOMES MODE ###################
    //######################################################################################
    
    else if (mode == "CompareGenomes"){
        std::shared_ptr<toyAgent> agent0 = std::make_shared<toyAgent>(&genomes[0]);
		agent0->RegulatoryPhenotype(toy);
		agent0->PatternFormation();
        std::shared_ptr<toyAgent> agent1 = std::make_shared<toyAgent>(&genomes[1]);
		agent1->RegulatoryPhenotype(toy);
		agent1->PatternFormation();

        if (agent1->PhenotypicPattern == agent0->PhenotypicPattern){
           cout << "Same pattern\n";
        }
        else{
            agent0->Print();
            agent1->Print();
        }

        return 0;
    }
	
	//######################################################################################
    //################################ STANDARD  MODE ######################################
    //######################################################################################

	/* ############## Setup ############## */
	
	Population* P;	

	/* ############## Initialisation ############## */
	printf("\033[93m### Initialisation ###\033[0m\n");
	P = new Population();
	//if(backup_reboot != "")	P->ContinuePopulationFromBackup();
	
	P->InitialisePopulation();
	
	printf("Initialisation completed...\n\n");

	

	/* ############## Simulation ############## */

	printf("\033[93m### Simulation ###\033[0m\n");
	for(Time=TimeZero; Time<SimTime+1; Time++){	//We do one extra step, because output is generated at the beginning of a step, such that time=0 is the field as it is initialised.
		

		
		if(Time%TimeSaveGrid==0)	P->PrintFieldToFile();
		if(Time%TimePruneFossils==0 && Time!=0)	P->PruneFossilRecord();
		if(Time%TimeOutputFossils==0)	P->Fossils->ExhibitFossils();	//Main next-state function, updating the population.

		if (Time%TimeTerminalOutput==0){
			P->Stats();
			cout << "Time: " << Time <<" Mean fitness: " << P->MeanFitness << " Max fitness: " << P->MaxFitness << " Size: " << P->size << endl;
		} 
		if (P->size==0) {
			printf("The Grim Reaper has taken over this simulation\n");
			break;
		}
		P->UpdatePopulation();		//Maybe move to main?????
		
	}

	/* ############## Final backups & ancestry ############## */

	Time--;
	
	P->PrintFieldToFile();
	if(SimTime%TimePruneFossils!=0)	P->PruneFossilRecord();
	if(SimTime%TimeOutputFossils!=0)	P->Fossils->ExhibitFossils();
	printf("Simulation completed...\n\n");

	/* ############## End ############## */

	printf("\033[93m### End ###\033[0m\n");
	delete P;
	P = NULL;
	printf("Prokaryotes completed...\n\n");

}




void Setup(int argc, char** argv) {

	string ReadOut, ReadOutN, command, name;
	bool project_name_found = false;
	bool initial_seed_set = false;
	


	for(int i=1;i<argc;i++)	//Loop through input arguments.
	{
		ReadOut = (char*) argv[i];	//There does not seem to be a quicker way to compare the input arguments with a string.

		//Let user define initial_seed. Otherwise defaults to time(0).
		if(ReadOut=="-s" && (i+1)!=argc)
		{
			initial_seed = atoi(argv[i+1]);
			initial_seed_set = true;
			printf("Seed = %i\n", initial_seed);
			i++;
			continue;
		}

		//Let user define subdirectory for the project
		else if(ReadOut=="-p" && (i+1)!=argc)
		{   
            name = argv[i+1];
			folder += name;
            OutputFile += name;
            OutputFile+=".txt";
			project_name_found = true;
			i++;
			continue;
		}

		else if(ReadOut=="-b" && (i+1)!=argc)
		{
			backup_reboot = argv[i+1];
			printf("Backup-file input: %s\n", backup_reboot.c_str());
			i++;
			continue;
		}

		else if(ReadOut=="-a" && (i+1)!=argc)
		{
			anctrace_reboot = argv[i+1];
			printf("Anctrace input: %s\n", anctrace_reboot.c_str());
			i++;
			continue;
		}

		else if(ReadOut=="-t" && (i+1)!=argc)
		{
			SimTime = atoi(argv[i+1]);
			i++;
			continue;
		}
		else if(ReadOut=="-tp" && (i+1)!=argc)
		{
			Target = atoi(argv[i+1]);
			i++;
			continue;
		}
		else if(ReadOut=="-ss" && (i+1)!=argc)
		{
			selection_strength = stod(argv[i+1]);
			printf("Selection parameter = %f\n", selection_strength);
			i++;
			continue;
		}
		else if(ReadOut=="--chrm")
		{
			chrom = true;
			continue;
		}
		else if(ReadOut=="--mix")
		{
			well_mixing = true;
			continue;
		}
		else if(ReadOut=="--sample")
		{
            mode = "GPsample";
			printf("Mode: genotype sample\n");
			continue;
		}
        else if(ReadOut=="--q")
		{
            mode = "quickSim";
			printf("Mode: quick simulations\n");
			continue;
		}
		else if(ReadOut=="-N" && (i+1)!=argc)
		{
			N = atoi(argv[i+1]);
			printf("N: %d\n", N);
			i++;
			continue;
		}
		else if(ReadOut=="--pg")
		{	
            mode = "printGenome";
            
			PrintGenome = true;
			printf("Mode: phenotype from genome\n");
			continue;
		}
		else if(ReadOut=="-g" && (i+1)!=argc)
		{
			gnm = argv[i+1];
			genome_init = std::string(argv[i+1]);
			printf("Genome: %s\n", genome_init.c_str());
			//printf("Genome: %s\n", anctrace_reboot.c_str());
			i++;
			continue;
		}
		else if(ReadOut=="-gs")
		{   
            mode = "CompareGenomes";
			int j=1;
			while (i+j != argc)
			{
				ReadOutN = (char*) argv[i+j];
				if (ReadOutN.substr(0,1)=="-" || j == 10000+1){
					i += j-1;	//The main for-loop will also add 1 to i.
					j=1;
					break;
				}
				else	//Haven't gotten to the next different command-line argument (or max. number of input files), so keep interpreting as another genome file.
				{
					genomes[j-1] = argv[i+j];
					j++;
				}
			}
			if (i+j == argc)
			{
				i += j-1;
				break;
			}
		}
		else	//Print usage/help.
		{
			printf("\n\033[93m### Patterns --- usage ###\033[0m\nArgument options:\n\t-p [project title]\n\t-s [seed]\n\t-g [initial genome]\n\t-b [backup file]\n\t-a [ancestor file]\n\t-t [max. time]\n\t-tp [target pattern]\n\t-ss [selection strength]\n\t--chrm [turn on chromosomal mutations]\n\t--mix [well mixing]\n\t--q [20 quick simulations]\nOther modes:\n\t--sample [gp sample]\n\t  -N [sample size]\n\t--pg [print genome data]\n\t  -g [genome or genome file] -gs [genomes for comparison]\n\n");
			exit(1);
		}
		
	}

	if (!initial_seed_set)	printf("Seed = %li\n", time(0));

	if (!project_name_found && mode=="standard")	folder += "You_are_lazy";	

	if (mode=="GPsample")	folder += "GP_Samples";	
	
	if (mode=="quickSim"){
		
		printf("NC: %d\n", NC);
		printf("NR: %d\n", NR);
		printf("Simulation time: %d\n", SimTime);
		printf("Mixing: %s\n", (well_mixing)? "yes" : "no");
		if (!well_mixing) printf("Diffusion rate: %f\n", diffusion_rate);
		printf("Initial genome size: %d\n", genome_length);
		printf("Mutation rates (point, delete, copy, cut): %f, %f, %f, %f\n", point_mut, del_mut, copy_mut, cut_mut);
		printf("Chromosomal mutations: %s\n", (chrom)? "yes" : "no");
		printf("Target pattern: %d\n", Target);
		
	}
	
	if (mode=="standard"){

		command = "mkdir -p " + folder;
		int r = system(command.c_str());
		printf("Folder = %s\n", folder.c_str());

		//Automatically set up a subdirectory for snapshots of the grid (not images but raw data).
		command = "mkdir -p " + folder + "/snapsamples";
		r += system(command.c_str());
		command = "mkdir -p " + folder + "/ancestors";
		r += system(command.c_str());
		if (r!=0) printf("Unable to create result directories\n");
		else printf("Result directories created succesfully...\n");
		printf("NC: %d\n", NC);
		printf("NR: %d\n", NR);
		printf("Simulation time: %d\n", SimTime);
		printf("Mixing: %s\n", (well_mixing)? "yes" : "no");
		if (!well_mixing) printf("Diffusion rate: %f\n", diffusion_rate);
		printf("Initial genome size: %d\n", genome_length);
		printf("Mutation rates (point, delete, copy, cut): %f, %f, %f, %f\n", point_mut, del_mut, copy_mut, cut_mut);
		printf("Chromosomal mutations: %s\n", (chrom)? "yes" : "no");
		printf("Target pattern: %d\n", Target);
		// for (const auto& row : Patterns[Target]) {
        //     for (int val : row) {
        //         cout << val << " ";
        //     }
        //     cout << "\n";
        // }
        // cout << "\n";
	}
	

}



