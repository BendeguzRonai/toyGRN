/*
* Code partly provided by Pablo Gonzalez Tamames
* https://github.com/tamames/toyLIFE_simulator
*/

#include "toyAgent.hh"
#include "Header.hh"

class DoubleComparator {
   public:
    bool operator()(const double lhs, const double rhs) const {
        // Check if the absolute difference is less than the threshold
        // Note: This is the comparison for ordering, not for equality
        // The map still requires a strict weak ordering
        double thr = 0.01;
        return lhs < rhs && std::abs(lhs - rhs) >= thr;
    }
};


toyAgent::toyAgent(string* genotype) {
    
    this->genotype = genotype;
    this->prots = {};
    this->dims = {};
    this->mets = {};
    this->vec_genotype = toy.genotype_str_to_vec(*genotype);
    this->g = vec_genotype.size();
    this->state_space = int_pow(2, g);
    this->logic_function= vector(3, std::vector<int>(state_space, 0));
    this->morphogen_outleak = std::vector<int>(state_space, 0);

    for (int p = 0; p < this->vec_genotype.size(); ++p) {

        
        int prom = this->vec_genotype[p].first;
        int prot =toy.prot_gen[this->vec_genotype[p].second]; 
        if (vocal) std::cout << "gene " << p << ": "<< prom << " " << prot  << "\n";
    }
}
toyAgent::~toyAgent() {


}


int toyAgent::promoter_expression(const ToyPlugin& toy) {
    // OUTPUT
    // PROMOTER_FIGHT GIVES THE OUTPUT STATE FROM A GIVEN SET OF OBJECTS
    // (PROTEINS (INT) AND DIMERS (BOOST::TUPLE(INT,INT,INT)) ATTENTION: IN THIS
    // VERSION OF PROMOTER_EXPRESSION, WE NEED NOT CHECK IF THE COMPETING
    // PROTEINS ARE THE SAME OR NOT AS THEY ARE ALREADY SAVED IN DIFFERENT SLOTS
    // IN THE OBJECTS MAP
    // VEC GENOTYPE
    // auto vgenotype = toy.genotype_str_to_vec(
    //     *genotype);  // transforms binary string into list of promoters and genes
    // LOOK UP THE BINDING ENERGY OF EVERY PROTEIN TO THE PROM

    // std::cout << "Genes " << std::endl;
    // for (auto it=vgenotype.begin(); it!=vgenotype.end(); ++it)
    //   std::cout << (*it).first << "\t" << toy.prot_gen[(*it).second] <<
    //   std::endl;
    // std::getchar();
    std::vector<bool> promoter_active(vec_genotype.size(), 0);

    // LOOP OVER ALL PROMOTERS
    for (int p = 0; p < vec_genotype.size(); ++p) {

        
        int prom = vec_genotype[p].first;
        if (vocal) std::cout<< "prom: " << prom << " poly: " << toy.polymerase[prom] << endl;
        //printMap(prots);
        if (vocal) std::cout << "gene " << p << ": "<< prom << " " << vec_genotype[p].second << " " << toy.prot_gen[vec_genotype[p].second]  << "\n";
        std::pair<std::map<Prot, double>, std::map<Dim, double>>
            binding_energies;
        for (mapa_prot::const_iterator it_prot = prots.begin();
             it_prot != prots.end(); ++it_prot){  // for every protein in PROTS
            
            //std::cout << "prot: "<< it_prot->first << endl;

            if (d_less(toy.prot_prom_energies[prom][it_prot->first],
                       0.0) &&  // if the protein can bind to PROM AND
                d_less(toy.prot_prom_energies[prom][it_prot->first] +
                           toy.prot_energies[it_prot->first] -
                           toy.polymerase[prom],
                       0.0)){  // it can bind it more strongly than POLY
                binding_energies.first[it_prot->first] =
                    toy.prot_prom_energies[prom][it_prot->first] +
                    toy.prot_energies[it_prot->first] -
                    toy.polymerase[prom];  // we subtract the binding energy of
                                           // the polymerase; if it is greater
                                           // than zero, binding energy is zero.
                    //std::cout << "bind E: " << toy.prot_prom_energies[prom][it_prot->first] << "prot E: " << toy.prot_energies[it_prot->first] << endl;
                }
        }
        for (mapa_dim::const_iterator it_dim = dims.begin();
             it_dim != dims.end(); ++it_dim) {
            if (vocal) std::cout<<"dim "<< it_dim->first.id<<endl;
            if (vocal) std::cout << "Prom: " << prom << " dimer binding E: " <<toy.dim_prom_energies[prom][it_dim->first.id] << endl;
            if (vocal) std::cout << toy.dim_prom_energies[prom][it_dim->first.id] +
                           toy.dim_bond_energy[it_dim->first.id] +
                           toy.prot_energies[it_dim->first.p1] +
                           toy.prot_energies[it_dim->first.p2] -
                           toy.polymerase[prom] << endl;
            if (d_less(toy.dim_prom_energies[prom][it_dim->first.id],
                       0.0) &&  // if the dimer can bind to PROM AND
                d_less(toy.dim_prom_energies[prom][it_dim->first.id] +
                           toy.dim_bond_energy[it_dim->first.id] +
                           toy.prot_energies[it_dim->first.p1] +
                           toy.prot_energies[it_dim->first.p2] -
                           toy.polymerase[prom],
                       0.0))  // it can bind it more strongly than POLY
                binding_energies.second[it_dim->first] =
                    toy.dim_prom_energies[prom][it_dim->first.id] +
                    toy.dim_bond_energy[it_dim->first.id] +
                    toy.prot_energies[it_dim->first.p1] +
                    toy.prot_energies[it_dim->first.p2] -
                    toy.polymerase[prom];  // we subtract the binding energy of
                                           // the polymerase; if it is greater
                                           // than zero, binding energy is zero
             } // for every dimer in DIMS
            
        // WE LOOK FOR THE MINIMUM
        while (1) {
            double min = 0.0;
            bool i_min = 0;  //"0" means the minimum is a protein; "1" means it
                             // is a dimer
            std::map<Prot, double>::const_iterator prot_min;
            std::map<Dim, double>::const_iterator dim_min;
            for (std::map<Prot, double>::const_iterator it1 =
                     binding_energies.first.begin();
                 it1 != binding_energies.first.end();
                 ++it1) {  // loop over all proteins, we save its id (i_min) and
                           // its energy (min)
                if (d_less(it1->second, min)) {
                    min = it1->second;
                    prot_min = it1;
                }  // if binding<min
            }  // for i
            for (std::map<Dim, double>::const_iterator it2 =
                     binding_energies.second.begin();
                 it2 != binding_energies.second.end();
                 ++it2) {  // loop over all dimers, we save its id (i_min) and
                           // its energy (min)
                if (d_less(it2->second, min)) {
                    min = it2->second;
                    i_min = 1;
                    dim_min = it2;
                }  // if binding<min
            }  // for i
            //std::cout << prom << " =prom | min=  " << min << endl;
            //printMap(prots);
            if (d_equal(min,
                        0.0)) {  // If the minimal is 0, no protein or dimer
                                 // binds, and it all depends on the polymerase
                //std::cout << "nth binded to prom " << prom << " " << d_less(toy.polymerase[prom], 0.0) <<"\n";
                promoter_active[p] = d_less(toy.polymerase[prom], 0.0);
                
                break;  // we exit the loop
            }

            // We check if there is more than one minimum, and eliminate the
            // repeated ones
            int count_mins = 0;
            if (!i_min) {  // if the minimum is a protein
                for (std::map<Prot, double>::const_iterator it1 =
                        std::next(prot_min);
                     it1 != binding_energies.first.end();
                     ++it1)  // we start the loop on the object following the
                             // minimal (the others have greater binding
                             // energies)
                    if (d_equal(it1->second,
                                min))  // if the binding energies are the same
                        count_mins++;  // we sum to the counter
                for (std::map<Dim, double>::const_iterator it2 =
                         binding_energies.second.begin();
                     it2 != binding_energies.second.end();
                     ++it2)  // we compare the minimum protein with the dimers
                    if (d_equal(it2->second, min))
                        count_mins++;
            }  // if minimum is a protein
            else {  // minimum is a dimer (then we don't need to check
                    // repetitions in proteins, because we examined the proteins
                    // first)
                for (std::map<Dim, double>::const_iterator it2 =
                         std::next(dim_min);
                     it2 != binding_energies.second.end();
                     ++it2)  // we start the loop on the object following the
                             // minimal (the others have greater binding
                             // energies)
                    if (d_equal(it2->second,
                                min))  // if the binding energies are the same
                        count_mins++;  // we sum to the counter
            }  // minimum is a dimer

            if (count_mins) {  // if there are repeated energies
                for (auto it1 =binding_energies.first.begin(); it1 != binding_energies.first.end();){
		  if (d_equal(it1->second,
			      min))  // if the binding energies are the same
		    it1=binding_energies.first.erase(it1);
		  else
		    ++it1;// we sum to the counter
		}
                for (auto it2 =binding_energies.second.begin(); it2 != binding_energies.second.end();){
		  if (d_equal(it2->second, min))
		    it2=binding_energies.second.erase(it2);
		  else
		    ++it2;
		}
	      
            } else {  // if there is only one minimum, we check if PROM is
                      // activated
                if (!i_min)
                    promoter_active[p] =
                        toy.prot_pol[prom]
                                    [prot_min->first];  // return the minimum
                else
                    promoter_active[p] =
                        toy.dim_pol[prom]
                                   [dim_min->first.id];  // return the minimum
                break;                                   // we exit the loop
            }
        }  // while TRUE
    }  // for promoters

    // EXPRESSION
    //printMap(prots);
    prots.clear();  // now we know which promoters are expressed, all objects
                    // disappear
    dims.clear();
    for (int i = 0; i < promoter_active.size(); ++i) {
        if (promoter_active[i]) {
            prots[toy.prot_gen[vec_genotype[i].second]]++;  // we add the proteins
                                                         // that are expressed
        }
    }
    

    // OUTPUT
    int translated = 0;  // as every gene has a place in the genotype, we can express every
            // state as an integer: (1,1,0) = 2^2+2^1 = 6
    for (int i = 0; i < promoter_active.size(); ++i){
    //std::cout << vec_genotype[promoter_active[i]].second << endl;
    if (vocal) std::cout  << "prom " << i << " is " <<  promoter_active[i] << " " <<promoter_active[i] * int_pow(2, i) << "\n";
        if (toy.prot_gen[vec_genotype[i].second]!=0) 
            translated +=  promoter_active[i] * int_pow(2, i);
    }
    return translated;
}

void toyAgent::reacting(const ToyPlugin& toy) {
    // GIVEN A SET OF OBJECTS, IT TELLS US WHAT BINDS WHAT AND WHAT BREAKS
    // WHAT
    std::map<double, mapa_owm, DoubleComparator> binding_energies;  // creates a map of energies, so that complexes
                           // with the same final energy are classified in
                           // the same group. Each energy has a map of OWM,
                           // which has three members: prot, dim and met (if
                           // there is no Prot in the complex, it will be
                           // -1; if there is no Dim, it will be (-1,0,0),
                           // and if there is no Met it will be the empty
                           // string, "") and [VALUE] an int counting how
                           // many possible molecules can be formed of this
                           // type
    // First we check if any Dimer forms
    for (mapa_prot::const_iterator it1 = prots.begin(); it1 != prots.end();
         ++it1)  // we only want to explore the combinations of two proteins
        for (mapa_prot::const_iterator it2 = it1; it2 != prots.end();
             ++it2)  // remember that proteins are sorted (and that they can
                     // form dimers with themselves)
            if (it1 != it2 ||
                (it1 == it2 &&
                 it1->second >
                     1)) {  // Prots will only dimerize with themselves if
                            // there is more than one copy
                std::pair<int, int> perim =
                    std::make_pair(toy.prot_perimeters[it1->first],
                                   toy.prot_perimeters[it2->first]);
                std::map<std::pair<int, int>, int>::const_iterator it_dim =
                    toy.dim_perim.find(perim);
                if (it_dim !=
                    toy.dim_perim.end()) {  // if this pair forms a dimer,
                                            // we form the tuple
                    double b_energy =
                        toy.dim_bond_energy[it_dim->second] +
                        toy.prot_energies[it1->first] +
                        toy.prot_energies[it2->first];  // energy of the
                                                        // complex
                    if (vocal) std::cout << "proteins: " << it1->first <<" & " << it2->first << " E: " << b_energy << "\n";
                    Dim new_dim(it_dim->second, it1->first, it2->first);
                    OWM new_complex(-1, new_dim, Met());
                    binding_energies[b_energy][new_complex] =
                        int_min(it1->second, it2->second);
                }  // if a dimer can be formed
            }  // if it1!=it2

    // Then we check if any Protein binds any Met
    // for (mapa_prot::const_iterator it_prot = prots.begin();
    //      it_prot != prots.end(); ++it_prot)
    //     for (mapa_met::const_iterator it_met = mets.begin();
    //          it_met != mets.end(); ++it_met) {
    //         pairPM par1 = std::make_pair(it_met->first, it_prot->first);
    //         try {  // see if the prot binds the met
    //             // std::map<Met, std::map<Prot,double> >::const_iterator
    //             // iter_tuplita=prot_met.find(it_met->first);
    //             // std::map<Prot,double> tuplita=iter_tuplita->second;
    //             double e1 = toy.prot_met.at(par1);
    //             if (d_less(e1,
    //                        0.0)) {  // if the protein can bind the metabolite
    //                 double b_energy =
    //                     e1 + toy.prot_energies[it_prot->first];  // energy of
    //                                                              // the complex
    //                 OWM new_complex(it_prot->first, Dim(), it_met->first);
    //                 binding_energies[b_energy][new_complex] =
    //                     int_min(it_prot->second, it_met->second);
    //             }
    //         } catch (std::out_of_range) {
    //         }  // if it doesn't, don't do anything
    //     }

    // Then we check if any Dimer binds any Met
    // for (mapa_dim::const_iterator it_dim = dims.begin(); it_dim != dims.end();
    //      ++it_dim)
    //     for (mapa_met::const_iterator it_met = mets.begin();
    //          it_met != mets.end(); ++it_met) {
    //         pairDM par1 = std::make_pair(it_met->first, it_dim->first.id);
    //         try {  // look and see if the dimer binds to the met
    //             Dmet tupla = toy.dim_met.at(par1);
    //             if (d_less(tupla.eg,
    //                        0.0)) {  // if the dimer can bind the metabolite
    //                 double b_energy =
    //                     tupla.eg + toy.dim_bond_energy[it_dim->first.id] +
    //                     toy.prot_energies[it_dim->first.p1] +
    //                     toy.prot_energies[it_dim->first.p2];  // energy of
    //                                                           // the complex
    //                 OWM new_complex(-1, it_dim->first, it_met->first);
    //                 binding_energies[b_energy][new_complex] =
    //                     int_min(it_dim->second, it_met->second);
    //             }
    //         } catch (std::out_of_range) {
    //         }  // if it doesn't don't do anything
    //     }

    // // Finally, we check if any D+M breaks
    // for (mapa_owm::const_iterator it_dm = owns.begin(); it_dm != owns.end();
    //      ++it_dm) {
    //     if (it_dm->first.dim.id == -1)
    //         continue;  // this means it's a P+M complex
    //     for (mapa_prot::const_iterator it_prot = prots.begin();
    //          it_prot != prots.end(); ++it_prot) {
    //         pairDM par1 = std::make_pair(it_dm->first.met, it_dm->first.dim.id);
    //         Dmet tupla = toy.dim_met.at(par1);
    //         if (tupla.seq.size() != 16)
    //             continue;  // only check this if the D+M can actually be
    //                        // broken!
    //         std::pair<int, std::string> par2(it_prot->first, tupla.seq);
    //         try {  // look and see if the protein can break this dimer
    //             double b_energy = toy.prot_breaking.at(par2).first +
    //                               toy.prot_energies[it_prot->first];
    //             int c_min = toy.prot_breaking.at(par2).second;
    //             if (c_min <
    //                 2)  // the sequences should be ordered in such a way
    //                     // that the first half is always the first protein
    //                 b_energy +=
    //                     toy.prot_energies[it_dm->first.dim.p1];  // sum energy
    //                                                              // of first
    //                                                              // subunit
    //             // if ((c_min<2 && tupla.bind<3) ||//if the prot binds to
    //             // the first subunit when the dimer is not reversed
    //             //     (c_min>1 && tupla.bind>2))//or if it binds to the
    //             //     second subunit when the dimer is reversed
    //             //   b_energy += prot_energies[it_dm->first.dim.p1];//sum
    //             //   energy of first subunit
    //             else
    //                 b_energy +=
    //                     toy.prot_energies[it_dm->first.dim.p2];  // sum energy
    //                                                              // of second
    //                                                              // subunit
    //             double met_energy =
    //                 toy.dim_bond_energy[it_dm->first.dim.id] +
    //                 toy.prot_energies[it_dm->first.dim.p1] +
    //                 toy.prot_energies[it_dm->first.dim.p2];  // E of the
    //                                                          // existing complex
    //             if (d_less(b_energy, met_energy)) {  // then it could break
    //                 OWM new_complex(it_prot->first, it_dm->first.dim,
    //                                 it_dm->first.met);
    //                 binding_energies[b_energy][new_complex] =
    //                     int_min(it_prot->second, it_dm->second);
    //             }
    //         }  // try
    //         catch (std::out_of_range) {
    //         }  // it doesn't bind to this sequence
    //     }  // for prots
    // }  // for d+m complexes

    // NOW WE HAVE A MAP OF ENERGIES
    // LOOP OVER EACH ENERGY AND SEE WHAT HAPPENS
    for (auto it = binding_energies.begin(); it != binding_energies.end();
         ++it) {
        mapa_prot number_prots;  // will hold how many complexes each
                                 // protein is involved in (if it's greater
                                 // than prots[i], then no dimers is formed
        mapa_dim number_dims;    // will hold how many complexes each protein
                                 // is involved in (if it's greater than
                                 // prots[i], then no dimers is formed
        mapa_met number_mets;    // will hold how many complexes each protein
                                 // is involved in (if it's greater than
                                 // prots[i], then no dimers is formed
        mapa_owm number_dm;      // because std::map doesn't work well with
                                 // boost::tuple, this vector will serve as an
                                 // index for the elements in owns;
        for (mapa_owm::iterator it_tuple = it->second.begin(); it_tuple != it->second.end(); ++it_tuple){  // for each complex
            
            if (it_tuple->first.met.empty()) {  // protein-protein complex
                number_prots[it_tuple->first.dim.p1] += it_tuple->second;

                if (it_tuple->first.dim.p2!=it_tuple->first.dim.p1) //bug that didn't allow homodimers fixed?
                number_prots[it_tuple->first.dim.p2] += it_tuple->second;
            } 
            else if (it_tuple->first.prot != -1) {  // there are proteins and mets
                if (it_tuple->first.dim.empty()) {  // protein-met complex (no dimers)
                    number_prots[it_tuple->first.prot] += it_tuple->second;
                    number_mets[it_tuple->first.met] += it_tuple->second;
                } 
                else {  // protein breaking a d+m complex
                    number_prots[it_tuple->first.prot] += it_tuple->second;
                    number_dm[OWM(-1, it_tuple->first.dim,
                                  it_tuple->first.met)] += it_tuple->second;
                }
            }  // proteins and mets are involved
            else {  // dimer-met complex
                number_dims[it_tuple->first.dim] += it_tuple->second;
                number_mets[it_tuple->first.met] += it_tuple->second;
            }
        }
        // NOW, FOR EVERY POSSIBLE COMPLEX, WE CHECK IF THE UNITS ARE
        // INVOLVED IN TOO MANY COMPLEXES AND DISCARD THOSE
        for (mapa_owm::iterator it_tuple = it->second.begin();
             it_tuple != it->second.end();) {   // for each complex
            if (it_tuple->first.met.empty()) {  // protein-protein complex
                
                //printf("we re here, p1: %d, count: %d | p2: %d, count: %d\n", it_tuple->first.dim.p1, prots[it_tuple->first.dim.p1], it_tuple->first.dim.p2, prots[it_tuple->first.dim.p2]);
                //printf("count in dimersL p1: %d, p2: %d\n", number_prots[it_tuple->first.dim.p1], number_prots[it_tuple->first.dim.p2]);
                if (prots[it_tuple->first.dim.p1] <
                        number_prots[it_tuple->first.dim.p1] ||  // prot1 more
                                                                 // involved
                    prots[it_tuple->first.dim.p2] <
                        number_prots[it_tuple->first.dim.p2]) {  // prot2 more
                                                                 // involved
                    it_tuple =
                        it->second.erase(it_tuple);  // this dimer will not form

                    //printf("we re here22, %d\n", it_tuple->second);
                } else
                    ++it_tuple;
            } else if (it_tuple->first.prot !=
                       -1) {  // proteins and mets are involved
                if (it_tuple->first.dim.empty()) {  // protein-met complex
                    if (prots[it_tuple->first.prot] <
                            number_prots[it_tuple->first.prot] ||  // proteins
                                                                   // more
                                                                   // involved
                        mets[it_tuple->first.met] <
                            number_mets[it_tuple->first.met])  // met more
                                                               // involved
                        it_tuple = it->second.erase(
                            it_tuple);  // this P+M will not form
                    else
                        ++it_tuple;
                } else {  // protein breaking a d+m complex
                    OWM new_tuple(-1, it_tuple->first.dim, it_tuple->first.met);
                    if (prots[it_tuple->first.prot] <
                            number_prots[it_tuple->first
                                             .prot] ||  // protein more involved
                        owns[new_tuple] <
                            number_dm[new_tuple])  // d+m more involved
                        it_tuple = it->second.erase(
                            it_tuple);  // this D+M will not be broken
                    else
                        ++it_tuple;
                }
            }  // proteins and mets are involved
            else {  // dimer-met complex
                if (dims[it_tuple->first.dim] <
                        number_dims[it_tuple->first.dim] ||  // dimer more
                                                             // involved
                    mets[it_tuple->first.met] <
                        number_mets[it_tuple->first.met])  // met more involved
                    it_tuple =
                        it->second.erase(it_tuple);  // this D+M will not form
                else
                    ++it_tuple;
            }
        }

        // NOW WE FORM THE COMPLEXES THAT REMAIN
        for (mapa_owm::const_iterator it_tuple = it->second.begin();
             it_tuple != it->second.end(); ++it_tuple) {  // for each complex
            if (it_tuple->first.met.empty()) {  // protein-protein complex
                if (it_tuple->first.dim.p1!=it_tuple->first.dim.p2){
                dims[it_tuple->first.dim] += it_tuple->second;  // we add as many copies as possible
                prots[it_tuple->first.dim.p1] -= it_tuple->second;  // we eliminate the Prot copies
                prots[it_tuple->first.dim.p2] -= it_tuple->second;
                }
                else {
                    int dim_n = it_tuple->second/2;
                    dims[it_tuple->first.dim] += dim_n;  // we add as many copies as possible
                    prots[it_tuple->first.dim.p1] -= dim_n*2;  // we eliminate the Prot copies
                    
                }

            } else if (it_tuple->first.prot !=
                       -1) {  // there are proteins and mets
                if (it_tuple->first.dim.empty()) {  // protein-met complex
                    owns[it_tuple->first] += it_tuple->second;
                    prots[it_tuple->first.prot] -= it_tuple->second;
                    mets[it_tuple->first.met] -= it_tuple->second;
                } else {  // protein breaking a d+m complex
                    owns[OWM(-1, it_tuple->first.dim, it_tuple->first.met)] -=
                        it_tuple->second;
                    prots[it_tuple->first.prot] -= it_tuple->second;
                      // it->first?
                    // Now we add the rests of the mets to Met
                    pairDM par1 = std::make_pair(it_tuple->first.met,
                                                 it_tuple->first.dim.id);
                    Dmet tupla = toy.dim_met.at(par1);
                    std::string pmet1 = it_tuple->first.met.substr(
                        0, tupla.pos);  // we split the met in the two parts
                                        // defined by the Dimer when it binds
                    std::string pmet2 =
                        it_tuple->first.met.substr(pmet1.size());
                    mets[pmet1] += it_tuple->second;
                    mets[pmet2] += it_tuple->second;
                }
            }  // proteins and mets
            else {  // dimer-met complex (it should check this: if
                // (boost::get<0>(boost::get<1>(*it_tuple))!=-1)
                owns[it_tuple->first] += it_tuple->second;
                dims[it_tuple->first.dim] -= it_tuple->second;
                mets[it_tuple->first.met] -= it_tuple->second;
            }
        }  // for each complex
    }  // for all energies

    // NOW WE ELIMINATE THE ELEMENTS THAT HAVE BEEN USED UP
    for (auto it_prot = prots.begin(); it_prot != prots.end();)
        if (it_prot->second == 0)
            it_prot = prots.erase(
                it_prot);  // eliminate the proteins that have 0 copies
        else
            ++it_prot;

    for (auto it_dim = dims.begin(); it_dim != dims.end();) {
        if (it_dim->second == 0)
            it_dim = dims.erase(it_dim);  // eliminate dimers that have 0 copies
        else
            ++it_dim;
    }

    for (auto it_met = mets.begin(); it_met != mets.end();)
        if (it_met->second == 0)
            it_met = mets.erase(it_met);  // eliminate mets that have 0 copies
        else
            ++it_met;

    for (auto it_owm = owns.begin(); it_owm != owns.end();)
        if (it_owm->second == 0)
            it_owm = owns.erase(it_owm);
        else
            ++it_owm;

    return;
}


// REGULATORY PHENOTYPE (COMPLETE BOOLEAN FUNCTION in toy plugin)

void toyAgent::RegulatoryPhenotype(const ToyPlugin& toy) {
    // COMPUTES THE ATTRACTOR OF THE DYNAMICS OF A toyLIFE ORGANISM

    std::vector<bool> visited(state_space, 0);
    int morphogen = toy.prot_gen[vec_genotype[0].second];
    if (vocal) printf("morph: %d \n", morphogen );
    //for (int it=0; it < 3; it++) std::cout << it << " " << toy.prot_gen[vec_genotype[it].second] << " ";
    //std::cout << endl;

    for (int i = 0; i < state_space; ++i){ //checking leak in main loop is buggy, will transfer here
        prots.clear();
        dims.clear();
        std::vector<int> who = str_to_vec(dectobin(i, g));
        std::reverse(who.begin(), who.end()); 

        if (who[0]){ // so that there is no leakage fromn identical genes
            for (int p = 0; p < g; ++p){
                if (who[p]) prots[toy.prot_gen[vec_genotype[p].second]]++;
            }
            reacting(toy); 

            if (prots[morphogen]>0) morphogen_outleak[i] = 1;
        }
        
        //std::cout << "state: " << i << " leak: " << morphogen_outleak[i] << endl;

    }

    for (int leak_in = 0; leak_in < 3; ++leak_in){
        //std::cout << leak_in << " " << logic_function[leak_in].size() << endl;

        for (int i = 0; i < state_space; ++i) {
        
            prots.clear();
            dims.clear();

            std::vector<int> who = str_to_vec(dectobin(i, g));
            std::reverse(who.begin(), who.end()); 
            //std::cout << "who: " << dectobin(i, g) << endl;
            if (vocal) printf("state: %d\n", i);
            
            for (int p = 0; p < g; ++p){
                
                if (vocal) std::cout << "gene: " << p << " " << toy.prot_gen[vec_genotype[p].second] << " ";
                
                if (who[p]){
                    prots[toy.prot_gen[vec_genotype[p].second]]++;
                    if (vocal) std::cout << "on ";
                }
                else if (vocal) std::cout << "off "; 
            }
            if (vocal) std::cout << endl;

            // std::cout << i << ": " << dectobin(i, g) << "\n";
           
            if (leak_in!= 0) prots[morphogen] += leak_in;
            
            reacting(toy);

            if (vocal){
                std::cout << "prots\n";
                printMap(prots);
                std::cout << "dims\n";
                printMap(dims);
            }
            
            
            int nu_state = promoter_expression(toy);  // compute next state
            
            if (vocal) {
                std::cout << "prots\n";
                printMap(prots);
                std::cout << "dims\n";
                printMap(dims);
            }

            logic_function[leak_in][i] = nu_state;
            
            if (vocal) cout << "new state: " << nu_state << "\n---------" <<endl;
        }
    }
}
int toyAgent::LogicAssert(){

    int error = 0;
    if (logic_function[0][1]!=logic_function[1][0]) error+=1; // 1xp0 
    if (logic_function[1][1]!=logic_function[2][0]) error+=10; // 2xp0
    if (logic_function[1][2]!=logic_function[0][3]) error+=100; // p0 & p1
    if (logic_function[1][3]!=logic_function[2][2]) error+=1000; // 2xp0 & p1
    
    return error;
}

void toyAgent::PatternFormation(){
    
    Pattern.resize(ncells, std::vector<int>(tmax, 0));
    
    int i, t, state, morphogen_infl;

    //Middle cell expresses A (Protein 01)
    int middle = ncells/2;
    Pattern[middle][0] = 1;
   
    // /*DEBUG*/
    // vector<int>::iterator it;
    // printf("Logic function:\n");
    // for(it = logic_function.begin(); it != logic_function.end(); it++, i++) {
	// 	// Increment distance for each differing character
	// 	printf("%d ", logic_function[i]);
	// }
    // printf("\nPattern:\n");

    // PrintPattern();

    for (t=0; t<tmax-1; t++){
        //std::cout << t << endl;
        morphogen_infl = morphogen_outleak[Pattern[1][t]];
        
        Pattern[0][t+1] = logic_function[morphogen_infl][Pattern[0][t]];
    
        morphogen_infl = morphogen_outleak[Pattern[29][t]];
        Pattern[30][t+1] = logic_function[morphogen_infl][Pattern[30][t]];

        for(i=1; i<ncells-1; i++){
            morphogen_infl = morphogen_outleak[Pattern[i-1][t]]+morphogen_outleak[Pattern[i+1][t]];
            
            Pattern[i][t+1] = logic_function[morphogen_infl][Pattern[i][t]];
            //if (t<2) std::cout << "index: " << i << " neigh state: " << Pattern[i-1][t] << ", " << Pattern[i+1][t] << " infl: " << morphogen_infl << " new state: " << Pattern[i][t+1] << endl;
        }
    }

    if (this->genotype->length()>59){
        PhenotypicPattern.resize(ncells, std::vector<int>(tmax, 0));
        
        for (t=0; t<tmax; t++){
            for (i=0; i<ncells; i++) PhenotypicPattern[i][t] = Pattern[i][t] %4;
        }
        
    }
    else {
      
        PhenotypicPattern = Pattern;
    }
}

double toyAgent::FitnessCalculation() {
    int h_dist;

    h_dist = CalculateHammingDistance();
    
    if (h_dist!=-1 && h_dist!=0) return pow((1.0-static_cast<double>(h_dist) / (ncells * tmax))*0.8, selection_strength);  
    else if (h_dist==0) return 1.;
    else return 0;

}

int toyAgent::CalculateHammingDistance() {

    bool nullPattern = true;
    for (int t = 0; t < tmax; ++t) {
        for (int k = 0; k < ncells; ++k) {
            if (Patterns[0][t][k] != PhenotypicPattern[k][t]) nullPattern = false;
            break;
        }
        if (!nullPattern) break;
    }
    if (nullPattern) return -1;
    int H=0;

    for (int t = 0; t < tmax; ++t) {
        for (int k = 0; k < ncells; ++k) {
            if (Patterns[Target][t][k] != PhenotypicPattern[k][t]) ++H;
        }
    }
    
    return H;
}

int toyAgent::FunctionToInt(){

    std::vector<int> func(8, -1);
    string s;
    int function;

    func[0] = this->logic_function[0][0];
    func[1] =this->logic_function[0][1] ;
    func[2] =this->logic_function[0][2];
    func[3] = this->logic_function[0][3] ;

    if (this->logic_function[1][0] != func[1]) printf("Error");//d
    func[4] = this->logic_function[1][1];
    if (this->logic_function[1][2] != func[3]) printf("Error"); // d
    func[5] = this->logic_function[1][3];

    if(this->logic_function[2][0] != func[4])  printf("Error"); //d
    func[6] = this->logic_function[2][1];
    if(this->logic_function[2][2] != func[5]) printf("Error"); //d
    func[7] = this->logic_function[2][3];

    if (this->morphogen_outleak[0] != 0) printf("Leak Error0\n"); 
    if(this->morphogen_outleak[1] != 1) printf("Leak Error1\n");  
    if(this->morphogen_outleak[2] != 0) {
        printf("Leak Error2\n");
        //this->Print();
    }
    int l = this->morphogen_outleak[3];//0 for dimer, 1 for leak

    for (int num : func) {
        s += std::to_string(num);
    }
    function = 100000*l + base4todec(s);
    return function;
}

void toyAgent::PrintPattern() {
    cout << "Regulatory Pattern\n";
    for (size_t t = 0; t < Pattern[0].size(); ++t) {  // Loop over columns (100 time steps)
        for (size_t i = 0; i < Pattern.size(); ++i) { // Loop over rows (31 cells)
            std::cout << Pattern[i][t] << " "; //" (" << PhenotypicPattern[i][t] << ") ";  // Print each element in the column
        }
        std::cout << std::endl;  // Move to the next line after printing each row (100 time steps)
    }

    cout << "Phenotypic Pattern\n";
    for (size_t t = 0; t < PhenotypicPattern[0].size(); ++t) {  // Loop over columns (100 time steps)
        for (size_t i = 0; i < PhenotypicPattern.size(); ++i) { // Loop over rows (31 cells)
            std::cout << PhenotypicPattern[i][t] << " "; //" (" << PhenotypicPattern[i][t] << ") ";  // Print each element in the column
        }
        std::cout << std::endl;  // Move to the next line after printing each row (100 time steps)
    }
    
    
}

void toyAgent::Print() {
    vector<int>::iterator it;
	int i = 0;
    std::cout << "Genotype: " << *genotype << "\n";
    
    
    printf("Logic function:\n");
    
    for (int leak_in=0; leak_in<3; leak_in++){
        i=0;
        printf("\nInfulx from neighbors: %d\n", leak_in);
        for(it = logic_function[leak_in].begin(); it != logic_function[leak_in].end(); it++, i++) 
        printf("%d ", logic_function[leak_in][i]);
    }
    printf("\nOutleakage:\n");
    
    i=0;
    for(it = morphogen_outleak.begin(); it != morphogen_outleak.end(); it++, i++) 
    printf("%d ", morphogen_outleak[i]);

    printf("\nGenes:\n");
    for (int p = 0; p < vec_genotype.size(); ++p) {
        printf("%d ",toy.prot_gen[vec_genotype[p].second]);
    }


    printf("\nPattern:\n");
    PrintPattern();
    
    // std::cout << "Prots: \n";
    // printMap(prots);
    // std::cout << "Dims: \n";
    // printMap(dims);
    // std::cout << "Mets: \n";
    // printMap(mets);
    // std::cout << "Owns: \n";
    // printMap(owns);
    
}
