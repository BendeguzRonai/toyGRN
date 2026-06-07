/*
* Based on Prokaryotes by Samuel von der Dunk
* https://github.com/samvonderdunk/Prokaryotes
*/
#include "Population.hh"


Population::Population()
{
	int i,j,n;
	p_nr_proks_=0;
	p_id_count_=0;

	Fossils = new FossilRecord();
	for(i=0;i<NR;i++) for(j=0;j<NC;j++)
	{
		PPSpace[i][j]=NULL;
	}
}

Population::~Population()
{
	int i,j,n;
	iterpps ips;

	for(i=0;i<NR;i++) for(j=0;j<NC;j++)
	{
		if((PPSpace[i][j])!=NULL)
		{
			if(PPSpace[i][j]->mutant)	Fossils->EraseFossil(PPSpace[i][j]->fossil_id);
			delete (PPSpace[i][j]);
			PPSpace[i][j]=NULL;
		}
	}
	
	delete Fossils;	//This will delete the rest of FossilList in the Fossils class internally.
	Fossils=NULL;
}

void Population::InitialisePopulation()
{
	Prokaryote* PP;
	Prokaryote* PP_Copy;

	//First create one Prokaryote.
	PP = new Prokaryote();
	PP->InitialiseProkaryote();
	PP->FitCalc();

	//Print its genome.
	//std::cout << "Initial genome = " << *(PP->G) << endl;	
	
	for(int row=0; row<NR; row++) for(int col=0; col<NC; col++){
	
		if (uniform() < 0.1)	//Initialise lower density
		{
			p_id_count_++;	//Make sure the first individual gets p_id_count_ of 1.
			
			//let's try to start from a diverse population
			if (DiverseInitialPop){ 
				delete PP;
				PP = NULL;
				PP = new Prokaryote();
				PP->InitialiseProkaryote();
				
			}
			
			PP_Copy=new Prokaryote();
			PP_Copy->Clone(PP, p_id_count_);
			PP_Copy->Ancestor = NULL;	//Null-pointer tells me the cell was initialised.
			PP_Copy->mutant = true; //to keep in fossil record
			PP_Copy->FitCalc();

			PPSpace[row][col] = PP_Copy;
			Fossils->BuryFossil(PPSpace[row][col]);
			//PP_Copy->PrintData(true);
			/*DEBUG*/
			//PPSpace[row][col]->PrintData(true);
		}
	}

	delete PP;
	PP = NULL;
}



/* ######################################################################## */
/*				UPDATE	UPDATE	UPDATE	UPDATE	UPDATE	UPDATE	UPDATE						*/
/* ######################################################################## */

void Population::UpdatePopulation()	//This is the main next-state function.
{
	int update_order[NR*NC];
	int u, i, j, e_block, s_block;
	double resource, diffusion_steps;
	Prokaryote* PP;
	std::list<Prokaryote*> MutantsList;

	for(u=0; u<NR*NC; u++) update_order[u]=u;
	random_shuffle(&update_order[0], &update_order[NR*NC]);		//Is also set by initial_seed through srand(initial_seed); see World.cc

	for(u=0; u<NR*NC; u++)		//Go through the field: birth, death (and later possibly diffusion).
	{
		i = update_order[u]/NC;	//Row index.
		j = update_order[u]%NC;	//Column index.
		
		if (PPSpace[i][j] != NULL and PPSpace[i][j]->time_of_appearance!=Time)	//Alive	 site.	Programme is written such that nothing happens at empty sites.
		{	
			if (uniform() < death_rate)
			{	
				/*DEBUG*/
				//printf("DEATH\nCell #%d, born at T=%d by #%d dies at T=%d\n", PPSpace[i][j]->fossil_id, PPSpace[i][j]->time_of_appearance, (PPSpace[i][j]->Ancestor!=NULL) ? PPSpace[i][j]->Ancestor->fossil_id : 0,Time);
				if (!PPSpace[i][j]->mutant) delete PPSpace[i][j];
				else PPSpace[i][j]->alive=false; //will be deleted when prunning fossil record	
				PPSpace[i][j] = NULL;
				continue;	
			}

			coords neigh = PickNeighbour(i, j);

			if (PPSpace[neigh.first][neigh.second]==NULL and uniform() < PPSpace[i][j]->fitness){//Only replicate on empty squares
				
				p_id_count_++;
				PP = new Prokaryote();
				
				PP->Clone(PPSpace[i][j], p_id_count_);
				PP->Mutate(); //if mutant, fitness becomes 0
				if (PP->mutant) MutantsList.push_back(PP); //fitnes & fossil record afterwards, so that cells don't reproduce when they are born
				
				PPSpace[neigh.first][neigh.second] = PP;

				/*DEBUG*/
				//printf("BIRTH\nCell #%d, born at T=%d by #%d gives birth to cell %d at T=%d\n", PPSpace[i][j]->fossil_id, PPSpace[i][j]->time_of_appearance, (PPSpace[i][j]->Ancestor!=NULL) ? PPSpace[i][j]->Ancestor->fossil_id : 0, PPSpace[neigh.first][neigh.second]->fossil_id,Time);
				//if (PP->mutant) printf("Child is MUTANT\n");
				if(PPSpace[i][j]->time_of_appearance==Time) printf("toa=T= %d\n", Time);
			}
		}
	}
	/*DEBUG*/
	//printf("exit%d\n", MutantsList.size());

	iterpps m = MutantsList.begin();
	while(m != MutantsList.end())
	{	
		(*m)->FitCalc();
		Fossils->BuryFossil((*m));
		m++;
	}

	//Margolus diffusion or mixing
	if (well_mixing)	WellMix();
	else{	
		diffusion_steps = diffusion_rate;
		while (diffusion_steps > 0.)
		{
			if (diffusion_steps >= 1.)
			{
				MargolusDiffusion();
				diffusion_steps -= 1.;
			}
			else
			{
				if (uniform() < diffusion_steps)	MargolusDiffusion();
				diffusion_steps -= 1.;	//Make sure that it is not slightly above zero (?)
			}
		}}


}

void Population::Stats(){
	int i, j;
	double H = 0.;
	double F = 0.;
	size = 0;
	LeastDistance = 1.;
	MaxFitness = 0. ; 
	for(int u=0; u<NR*NC; u++)	
	{
		i =u/NC;	//Row index.
		j = u%NC;	//Column index.
		if (PPSpace[i][j] != NULL) {
			H+=PPSpace[i][j]->relative_distance;
			F+=PPSpace[i][j]->fitness;
			size++;
			if (PPSpace[i][j]->relative_distance < LeastDistance) LeastDistance = PPSpace[i][j]->relative_distance;
			if (PPSpace[i][j]->fitness > MaxFitness) MaxFitness = PPSpace[i][j]->fitness;
		}
	}

	if (size != 0)	{
		MeanDistance = static_cast<double>(H) / size;
		MeanFitness = static_cast<double>(F) / size;
	}
	else {
		MeanDistance = 0;
		MeanFitness = 0;
	}
		
}

void Population::MargolusDiffusion()	//Based on Brem's function from Evolvabear_V9.0
{
	int kernel, i, j;
	Prokaryote* PPtemp;

	//Rotation of 2x2 squares divided up according to two different kernels (shift by one), and both run the entire field.
	//kernel=0 does not need to deal with wrapped boundaries if NR and NC are even, but kernel=1 does.
	for (kernel=0; kernel<2; kernel++)
	{
		for(i=kernel; i<NR-1+kernel; i+=2)	for(j=kernel; j<NC-1+kernel; j+=2)
		{
			if (uniform()<0.5)	//Counter-clockwise.
			{
				PPtemp = PPSpace[(i+1)%NR][j%NC];	//The %NR and %NC make sure that the grid is wrapped (e.g. the boundary row index "100" will return 100%100 = 0 if NR=100).
				PPSpace[(i+1)%NR][j%NC] = PPSpace[i%NR][j%NC];
				PPSpace[i%NR][j%NC] = PPSpace[i%NR][(j+1)%NC];
				PPSpace[i%NR][(j+1)%NC] = PPSpace[(i+1)%NR][(j+1)%NC];
				PPSpace[(i+1)%NR][(j+1)%NC] = PPtemp;
				PPtemp = NULL;
			}
			else								//Clockwise.
			{
				PPtemp = PPSpace[i%NR][(j+1)%NC];
				PPSpace[i%NR][(j+1)%NC] = PPSpace[i%NR][j%NC];
				PPSpace[i%NR][j%NC] = PPSpace[(i+1)%NR][j%NC];
				PPSpace[(i+1)%NR][j%NC] = PPSpace[(i+1)%NR][(j+1)%NC];
				PPSpace[(i+1)%NR][(j+1)%NC] = PPtemp;
				PPtemp = NULL;
			}
		}
	}
}

void Population::WellMix()
{
	int Space1D[NR*NC];
	int u, i, j;
	Prokaryote* SpaceMirror[NR][NC];

	for(u=0; u<NR*NC; u++)	Space1D[u]=u;
	random_shuffle(&Space1D[0], &Space1D[NR*NC]);
	for(u=0; u<NR*NC; u++)	SpaceMirror[u/NC][u%NC] = PPSpace[Space1D[u]/NC][Space1D[u]%NC];
	for(i=0; i<NR; i++) for(j=0; j<NC; j++)
	{
		PPSpace[i][j] = SpaceMirror[i][j];
		SpaceMirror[i][j] = NULL;	//Should not be really necessary.
	}
}

Population::coords Population::PickNeighbour(int i, int j)	//Pick random neighbour.
{
	int nrow = i, ncol = j, ni, nj, random_neighbour;

	while (nrow == i && ncol == j)	//Try again if you pick yourself.
	{
		random_neighbour = (int)(uniform()*9);
		ni = random_neighbour/replication_neighbourhood;
		nj = random_neighbour%replication_neighbourhood;

		//Wrap grid boundaries
		nrow = i+ni-1;
		if(nrow < 0)	nrow += NR;
		else if(nrow >= NR)	nrow -= NR;

		ncol = j+nj-1;
		if(ncol < 0)	ncol += NC;
		else if(ncol >= NC)	ncol -= NC;
	}
	return std::make_pair(nrow, ncol);
}

int Population::NeighbourhoodDensity(int i, int j)
{
	int ii, jj, nrow, ncol, density=0;

	for (ii=i-1; ii<=i+1; ii++) for (jj=j-1; jj<=j+1; jj++)
	{
		if (ii == i && jj == j)	continue;

		if (ii < 0)	nrow = ii + NR;	//-1 becomes -1+100=99, i.e. the last index of a row with 100 sites (0-99).
		else if (ii >= NR)	nrow = ii - NR;
		else	nrow = ii;
		if (jj < 0)	ncol = jj + NC;
		else if (jj >= NC)	ncol = jj - NC;
		else	ncol = jj;

		if (PPSpace[nrow][ncol] != NULL)	density++;
	}

	return density;
}

/* ######################################################################## */
/*				FOSSILS	FOSSILS	FOSSILS	FOSSILS	FOSSILS	FOSSILS	FOSSILS						*/
/* ######################################################################## */

void Population::PrintFieldToFile()
{
	FILE* f;
	char OutputFile[800];
	sprintf(OutputFile, "%s/snapsamples/field%08d.txt", folder.c_str(), Time);
	f=fopen(OutputFile, "w");
	if (f == NULL){	printf("Failed to open file for writing the field.\n");	}

	fprintf(f, "i\tj\tID\tAncID\tToA\tfit\tGnm\n");
	for (int i=0; i<NR; i++) for(int j=0; j<NC; j++) {	//Don't print row and col numbers to save memory, these can be extracted by secondary scripts.
		if(PPSpace[i][j]==NULL){
		 	fprintf(f, "%d\t%d\tNULL\n", i, j);
		}
		else{		//Print relevant variables to make snapshots.
			fprintf(f, "%d\t%d\t%llu\t%llu\t%d\t%f\t%s\n", i, j, PPSpace[i][j]->fossil_id, (PPSpace[i][j]->Ancestor!=NULL) ? PPSpace[i][j]->Ancestor->fossil_id : 0, PPSpace[i][j]->time_of_appearance, PPSpace[i][j]->fitness, PPSpace[i][j]->G->c_str());
			
		}
	}
	fclose(f);
}


void Population::PruneFossilRecord()
{
	std::list<unsigned long long> AllFossilIDs;
	int fossil_record_size;
	for(int i=0; i<NR; i++)	for(int j=0; j<NC; j++)
	{
		if(PPSpace[i][j]!=NULL)
		{
			// Last common ancestor of living individual located
			Prokaryote* lastCA = PPSpace[i][j]->Ancestor;
			
			while(lastCA != NULL)
			{	
				AllFossilIDs.push_back(lastCA->fossil_id);	// Added to list of agents we must keep
				lastCA = lastCA->Ancestor;
			}
		}
	}
	// Delete duplicates (e.g. Agent 9 in example asci will be located 4 times. Agent 11 two times, etc.)
	AllFossilIDs.sort();
	AllFossilIDs.unique();

	// Delete all in FossilList that are not in AllFossilIDs (unless they are still living):
	fossil_record_size = (*Fossils).FossilList.size();
	iterpps ip = (*Fossils).FossilList.begin();
	while(ip != (*Fossils).FossilList.end())
	{
		// Search if stored agent was also found by tracing back:
		unsigned long long fossilID = (*ip)->fossil_id;
		iterull findit = std::find(AllFossilIDs.begin(),AllFossilIDs.end(),fossilID);
		// If not, delete the fossil unless it is still alive or is still saved in the graveyard. If a prokaryote dies, the graveyard-flag remains for one ShowGeneralProgress() cycle at most, so that the fossil can be deleted at the next pruning step. If ShowGeneralProgress() precedes PruneFossilRecord(), this is issue is even avoided, because flags are already removed off dead prokaryotes.
		if(findit==AllFossilIDs.end() && !(*ip)->alive)
		{
			delete *ip;
			ip = Fossils->FossilList.erase(ip);
		}
		else
		{
			++ip;
		}
	}

	AllFossilIDs.clear();
	//std::cout << "ID count (" << p_id_count_ << ")\tFossil record (pruned from " << fossil_record_size << " to " << (*Fossils).FossilList.size() << ")" << endl;
}


/* ######################################################################## */
/*									READ BACKUP								*/
/* ######################################################################## */


/*
//Soon, priviliges will have to be read from the backup-file as well.
void Population::ContinuePopulationFromBackup()
{
	ReadBackupFile();
	if(anctrace_reboot != "")	ReadAncestorFile();	//Currently, the fossil_ids are missing from the backup-file so it is impossible to link the fossils to live prokaryotes. But in the new version this will be possible.
	else	//Otherwise we will reset the fossil_ids to prevent these getting out of hand.
	{
		p_id_count_ = 0;
		for (int i=0; i<NR; i++)	for(int j=0; j<NC; j++)
		{
			if (PPSpace[i][j]!=NULL)
			{
				p_id_count_++;
				PPSpace[i][j]->fossil_id = p_id_count_;	//Other things such as time_of_appearance and Ancestor are set to zero by the EmptyProkaryote function.
			}
		}
	}

	OutputBackup();
	PruneFossilRecord();

	if (environmental_noise)	NoiseEnvironment();
}

void Population::ReadBackupFile()
{
	ifstream infile(backup_reboot.c_str());
	string line, data;
	char* data_element;
	string::iterator sit;
	Genome::iter it;
	int i, read_integer = 0, index, begin_data, end_data, counter, success, stage, pfork, panti_ori, temp_is_mutant, temp_priv, nr=NR, nc=NC, init_seed, read_header = 0;	//For now, set nr and nc for backup-file to NR and NC parameters (i.e. if backup-file does not contain header; for old backups). May be removed later.
	bool is_mutant, priv;
	unsigned long long prok_id, sdraws;
	double deficit;
	Prokaryote* PP;
	Gene* gene;

	if (!infile.is_open())
	{
		printf("Backup-file could not be opened.\n");
		exit(1);
	}

	printf("Reading backup from file: %s\n", backup_reboot.c_str());
	int count_lines = 0;
	while(getline(infile,line))
	{

		//Read "Header" of backup file.
		if (line=="### Header ###")
		{
			read_header = 1;
			continue;
		}
		else if (line=="#### Main ####")
		{
			read_header = -1;
			continue;
		}
		else if (read_header==1)
		{
			data_element = (char*)line.c_str();
			success = sscanf(data_element, "NR:%d\tNC:%d", &nr, &nc);
			if(success != 2)
			{
				cerr << "Could not read NR and NC from backup-file.\n" << endl;
				exit(1);
			}
			else	cout << "nr=" << nr << ", nc=" << nc << endl;
			read_header = 2;
			continue;
		}
		else if (read_header==2)
		{
			data_element = (char*)line.c_str();
			success = sscanf(data_element, "Initial seed:%d\tSeed draws:%llu", &init_seed, &sdraws);
			if (success != 2)
			{
				cerr << "Could not read seed status from backup-file.\n" << endl;
				exit(1);
			}
			else
			{
				cout << "Simulating " << sdraws << " random draws, initial seed=" << init_seed << endl;
				for(i=0; (unsigned long long)i<sdraws; i++){
					uniform();
				}
			}
			read_header = 3;	//Not doing anything with this yet.
			continue;
		}

		//Start reading "Main" data from backup file.
		//Make sure that field dimensions of backup and current simulation are related to one another.
		if (count_lines/nc >= NR || count_lines%nc >= NC)
		{
			count_lines++;
			continue;		//Sites that lie on [x,y] in the backup-sim, but where either x or y is not within the limits of the current sim; we simply lose these individuals.
		}
		if (count_lines%10000 == 0 && count_lines!=0)	printf("%d\n", count_lines);	//Print some progress.
		if(line == "0")	PPSpace[count_lines/nc][count_lines%nc] = NULL;	//Empty sites.
		else
		{
			//Start new individual.
			PP = new Prokaryote();
			PP->EmptyProkaryote();
			PP->G->BeadList = new list<Bead*>();
			PP->G->GeneTypes = new vector<int>();
			PP->G->GeneStates = new vector<int>();

			//Read BeadList. Do this first, because it first sets GeneStates and GeneTypes (randomised or based on other input files). We can reset these below.
			begin_data = line.find_first_of("(");
			end_data = line.find_last_of(")");
			data = line.substr(begin_data, end_data-begin_data+1);
			PP->G->ReadBeadsFromString(data);

			//Read GeneStates.
			index = line.find("]");
			data = line.substr(1, index-1);
			counter = 0;
			sit = data.begin();
			while(sit != data.end())
			{
				if(*sit == ',')
				{
					PP->G->GeneStates->push_back(read_integer);	//Save value just read.
					counter++;
					read_integer = 0;
				}
				else if(*sit != ' ')
				{
					read_integer *= 10;
					read_integer += (int)*sit - 48;
				}
				sit++;
			}
			PP->G->GeneStates->push_back(read_integer);
			read_integer = 0;

			//Read GeneTypes.
			begin_data = line.find("[", index);
			end_data = line.find("]", index+1);
			data = line.substr(begin_data+1, end_data-begin_data-1);
			counter = 0;
			sit = data.begin();
			while(sit != data.end())
			{
				if(*sit == ',')
				{
					PP->G->GeneTypes->at(counter) = read_integer;	//Save value just read.
					counter++;
					read_integer = 0;
				}
				else if(*sit != ' ')
				{
					read_integer *= 10;
					read_integer += (int)*sit - 48;
				}
				sit++;
			}
			PP->G->GeneTypes->at(counter) = read_integer;
			read_integer = 0;

			//Read prokaryote data.
			begin_data = line.find_last_of("[");
			end_data = line.find_last_of("]");
			data = line.substr(begin_data+1, end_data-begin_data-1);
			data_element = strtok((char*)data.c_str(),"\t");
			while(data_element != NULL)
			{
				success = sscanf(data_element, "%d %lf %d %d %llu %d %d", &stage, &deficit, &pfork, &panti_ori, &prok_id, &temp_is_mutant, &temp_priv);	//WARNING: in newer backup-files, the AncID has been inserted before temp_is_mutant and temp_priv, so BEWARE!!!
				if(success != 7)
				{
					cerr << "Could not find sufficient information for this prokaryote. Backup file potentially corrupt.\n" << endl;
					exit(1);
				}
				is_mutant = temp_is_mutant;	//Integers to bools (maybe not completely necessary..?)
				priv = temp_priv;

				data_element = strtok(NULL, "\t");
				PP->Stage = stage;
				PP->fitness_deficit = deficit;
				PP->G->pos_fork = pfork;
				PP->G->pos_anti_ori = panti_ori;
				PP->fossil_id = prok_id;
				PP->mutant = is_mutant;
				PP->priviliges = priv;
				if (PP->Stage>=2)	PP->time_replicated = 1;	//NOTE: we're not doing anything with time_replicated currently, so it's fine to set everything to 1 (allowing mitosis).
				if (prok_id > p_id_count_) p_id_count_ = prok_id;
			}

			// Read expression of individual genes. This will work in the new version, where the nr of genes actually corresponds to the length of the gene expression data.
			begin_data = line.find("{");
			end_data = line.find("}");
			data = line.substr(begin_data+1, end_data-begin_data-2);
			sit = data.begin();
			it = PP->G->BeadList->begin();
			while (sit != data.end())
			{
				if(*sit != ' ')
				{
					while(!PP->G->IsGene(*it))
					{
						it++;	//Go through beads until you hit the next gene.
					}

					gene = dynamic_cast<Gene*>(*it);
					gene->expression = (int)*sit - 48;
				}
				else	it++;
				sit++;
			}

			PP->G->SetClaimVectors();

			PPSpace[count_lines/nc][count_lines%nc] = PP;
			if (PP->mutant)	Fossils->BuryFossil(PP);
		}
		if(count_lines<generation_sample)
		{
			if(PPSpace[count_lines/nc][count_lines%nc] != NULL)	PPSpace[count_lines/nc][count_lines%nc]->saved_in_graveyard = true;
			OldGeneration[count_lines] = PPSpace[count_lines/nc][count_lines%nc];	//Put a subset of prokaryote pointers in the OldGeneration array.
		}
		count_lines++;
	}
	if (count_lines!= nr*nc)
	{
		cout << count_lines << " " << nr << " " << nc << endl;
		printf("Length of backup file conflicts with header (nr*nc).\n");
		exit(1);
	}
}

void Population::ReadAncestorFile()
{
	ifstream infile(anctrace_reboot.c_str());
	string line, data;
	int begin_data, end_data, TimeOA;
	unsigned long long ID, AncID;
	iterpps ip, ip2;
	Prokaryote* PP;

	if (!infile.is_open())
	{
		printf("Ancestor-file could not be opened.\n");
		exit(1);
	}

	printf("Reading ancestors from file: %s\n", anctrace_reboot.c_str());
	int count_alive = 0;
	int count_fossils = 0;
	int count_lines = 0;
	while(getline(infile,line))
	{
		if (count_lines%10000 == 0 && count_lines!=0)	printf("%d\n", count_lines);
		end_data = line.find("\t");
		data = line.substr(0,end_data);
		stringstream(data) >> ID;

		begin_data = end_data;
		end_data = line.find("\t",end_data+1);
		data = line.substr(begin_data, end_data-begin_data);
		stringstream(data) >> AncID;

		begin_data = end_data;
		end_data = line.find("\t",end_data+1);
		data = line.substr(begin_data, end_data-begin_data);
		stringstream(data) >> TimeOA;

		begin_data = end_data;
		end_data = line.size();
		data = line.substr(begin_data+1, end_data-begin_data);

		ip = Fossils->FossilList.begin();
		while (ip != Fossils->FossilList.end())
		{
			if ((*ip)->fossil_id == ID)	//Then we have found a live prokaryote in our ancestor file, because it will have to be added to the FossilRecord from the backup file. For these guys we only have to find its ancestor in the FossilList (all other data has been read from the backup-file).
			{
				count_alive++;
				(*ip)->time_of_appearance = TimeOA;
				if (AncID == 0)	(*ip)->Ancestor = NULL;
				else
				{
					ip2 = Fossils->FossilList.begin();
					while(ip2 != Fossils->FossilList.end())
					{
						if ((*ip2)->fossil_id == AncID)
						{
							(*ip)->Ancestor = *ip2;
							break;
						}
						ip2++;
					}
					if(ip2 == Fossils->FossilList.end())
					{
						printf("Error: ancestor %llu not found...exiting.\n", ID);
						exit(1);
					}
				}
				break;
			}
			ip++;
		}

		if (ip == Fossils->FossilList.end())	//We did not break out of the loop, so we have apparently not encountered this ID among the current list of fossils.
		{
			count_fossils++;
			PP = new Prokaryote();
			PP->EmptyProkaryote();
			PP->fossil_id = ID;
			PP->alive = false;
			PP->mutant = true;
			PP->time_of_appearance = TimeOA;
			if(AncID == 0)	PP->Ancestor = NULL;
			else	//We have to find the rightful parent of this creature. It should be in the list because we read in the fossil record starting with the oldest fossils; everything afterwards should have a parent present in the record.
			{
				ip2 = Fossils->FossilList.begin();
				while (ip2 != Fossils->FossilList.end())
				{
					if ((*ip2)->fossil_id == AncID)
					{
						PP->Ancestor = *ip2;
						break;
					}
					ip2++;
				}
				if(ip2 == Fossils->FossilList.end())
				{
					printf("Error: ancestor %llu not found...exiting.\n", ID);
					exit(1);
				}
			}

			//Set up its ghost genome.	It only has beads.
			PP->G->BeadList = new list<Bead*>();
			PP->G->GeneTypes = new vector<int>();
			PP->G->GeneStates = new vector<int>();
			PP->G->ReadBeadsFromString(data);
			PP->G->pos_anti_ori = PP->G->g_length;	//Otherwise the function PrintContent(NULL, false, true) (i.e. printing only the parental genome to a file) while think the parental genome is non-existent (pos_anti_or = 0).
			Fossils->BuryFossil(PP);
		}
		count_lines++;
	}
}

*/