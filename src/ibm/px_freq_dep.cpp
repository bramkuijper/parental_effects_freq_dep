// Gillespie simulation 
// for the evolution of nongenetic effects
// in a hawk-dove game
// (c) Bram Kuijper 2016

#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <cassert>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "bramauxiliary.h"

//#define NDEBUG

using namespace std;


// random number generator 
// see http://www.gnu.org/software/gsl/manual/html_node/Random-Number-Generation.html#Random-Number-Generation 
gsl_rng_type const * T; // gnu scientific library rng type
gsl_rng *r; // gnu scientific rng 

const size_t Npatch = 4000; 
const size_t Npp = 2;
unsigned long int NumGen = 1e20; // maximum number of generations
int sample = 20;
double d0 = 0.01; //some starting values
double h0 = 0.01;

int seed = -1; // the random seed
time_t total_time;  // track the total time the simulation is running
unsigned long int generation = 0; // track the current generation 

// output once every 10^6 generations
unsigned long int skip = 1e6;

double d = 0.1 // dispersal probability

double mh1 = 1.0; // mortality of hawks when in a patch with a dove
double mh2 = 2.0; // mortality of hawks when in a patch with another hawk
double md1 = 1.0; // mortality of dove when in a patch with a dove
double md0 = 2.0; // mortality of dove when in a patch with another dove

// mutation rates
double mu_h = 0.01;
double sdmu_h = 0.01;

struct Patch {
    Individual locals[Npp];
};


Patch MetaPop[Npatch];

// statistics for sampling
// with 0, 1 or 2 hawks 
size_t Npatches[3];

// next to counts store lists with the ids
// of the corresponding patches
size_t patch_ids[3][Npatch];

// give the outputfile a unique name
string filename("sim_px");
string filename_new(create_filename(filename));
ofstream DataFile(filename_new.c_str());  // output file 

// initialize the command line arguments to vary 
// the parameters
void init_arguments(int argc, char *argv[])
{
    d = atof(argv[1]);
    v = atof(argv[2]);
    c = atof(argv[3]);
}

// initialization function, runs at start of sim
// gives all the individuals some starting
// value
void init_pop()
{ 
    // start the time
    total_time = time(NULL);

    // obtain a seed from current nanosecond count
	seed = get_nanoseconds();

    // set up the random number generators
    // (from the gnu gsl library)
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, seed);

    // initialize patch state counters and set them to 0
    for (int n_hawk = 0; n_hawk <= 2; ++n_hawk)
    {
        Npatches[n_hawk] = 0;
    }

    size_t n_hawk;

    // initialize all the patches
    for (size_t i = 0; i < Npatch; ++i)
    {
        n_hawk = 0;

        //initialize breeders
        for (size_t j = 0; j < Npp; ++j)
        {
            if (gsl_rng_uniform(r) < 0.5)
            {
                // randomly assign initial phenotypes
                MetaPop[i].locals[j].hawk = true;
                ++n_hawk
            }
            else
            {
                // randomly assign initial phenotypes
                MetaPop[i].locals[j].hawk = false;
            }

            // initialize allelic values
            for (int allele = 0; allele < 2; ++allele)
            {
                MetaPop[i].locals[j].phh[allele] = h0;
                MetaPop[i].locals[j].pdh[allele] = h0;
            }
        }

        // increase the counter of patches 
        // containing n_hawk individuals
        //
        size_t Npatches_of_this_type = 
            Npatches[n_hawk];

        // append this patch to the corresponding stack of patch id's
        patch_ids[n_hawk][Npatches_of_this_type] = i;

        // update counter
        ++Npatches[n_hawk];
    }
}


double mut_h(double p)
{
    assert(p >= 0 && p <= 1);
    p+= gsl_rng_uniform(r) < mu_h ? gsl_ran_gaussian(r,sdmu_h) : 0;
    p = p < 0 ? 0 : p > 1 ? 1 : p;

    return(p);
}

// make a kid out of parent asexually
// add it either to local offspring pool
// or to pool of dispersers
void create_kid(Individual &mother, Individual &Kid)
{
    // if it is a hawk mother, 
    // we'll need to express phh
    // otherwise pdh

    // mother's phenotype is one
    double expr_h = mother.hawk ? 0.5 * (mother.phh[0] + mother.phh[1]) 
                                :
                                0.5 * (mother.pdh[0] + mother.pdh[1]);

    // offspring gets hawk or dove phenotype
    Kid.hawk = gsl_rng_uniform(r) < expr_h;

    // inherit alleles
    for (int allele_i = 0; allele_i < 2; ++allele_i)
    {
        Kid.phh[allele_i] = mut_h(mother.phh[allele_i]);
        Kid.pdh[allele_i] = mut_h(mother.pdh[allele_i]);
    }
}

// do the stats on the data
void write_data()
{
    // get variance and means
    double phh,pdh,d,vard,varphh,varpdh;
    double meanphh = 0;
    double meanpdh = 0;
    double ssphh = 0;
    double sspdh = 0;

    // loop through patches and individuals
    // and get stats on genotypes and patch freqs.
    for (size_t i = 0; i < Npatch; ++i)
    {
        for (size_t j = 0; j < Npp; ++j)
        {
            phh = 0.5 * (MetaPop[i].locals[j].phh[0] + MetaPop[i].locals[j].phh[1]);
            pdh = .5 * (MetaPop[i].locals[j].pdh[0] + MetaPop[i].locals[j].pdh[1]);

            meanphh+=phh;
            meanpdh+=pdh;
            ssphh+=phh*phh;
            sspdh+=pdh*pdh;
        }    
    }

    DataFile << generation << ";"; 

    // write down all the patch frequency combinations
    for (size_t n_hawk = 0; n_hawk <= 2; ++n_hawk)
    {
        DataFile << setprecision(5) 
            << (double)Npatches[n_hawk] / Npatch << ";";
    }

    meanphh /= Npatch * 2;
    meanpdh /= Npatch * 2;
    varphh = ssphh / (Npatch * 2) - meanphh*meanphh;
    varpdh = sspdh / (Npatch * 2) - meanpdh*meanpdh;

    DataFile 
        << setprecision(5) << meanphh << ";"
        << setprecision(5) << meanpdh << ";"
        << setprecision(5) << varphh << ";"
        << setprecision(5) << varpdh << ";" << endl;
}

// headers of the datafile
void write_data_headers()
{
    DataFile << "generation;";

    for (size_t n_hawk = 0; n_hawk <= 2; ++n_hawk)
    {
        DataFile << setprecision(5) << "f_" << (n_hawk + 1) << ";";
    }
    
    DataFile << "meanphh;meanpdh;varphh;varphd;" << endl;
}

// parameters at the end of the sim
void write_parameters()
{
    total_time = time(NULL) - total_time;
    DataFile << endl << endl 
                << "d;" << d << endl
                << "v;" << v << endl
                << "c;" << c << endl
                << "mu_h;" << mu_h << endl
                << "sdmu_h;" << sdmu_h << endl
                << "NumGen;" << NumGen << endl
                << "seed;" << seed << endl
                << "Npatches;" << Npatches << endl
                << "runtime;" << total_time << endl;
}


// mortality of an individual in one of the patches 
//
// size_t n_hawk: a patch containing n_hawk adapted individuals
// bool mortality_maladapted: true, maladapted guy dies; false, adapted guy dies;

void mortality(size_t patch_envt, size_t n_hawk, bool mortality_hawk)
{
    // sample a random patch that fulfills the conditions
    size_t random_patch = 
        gsl_rng_uniform_int(r, 
            Npatches[n_hawk]
        );

    // get the patch_id 
    size_t random_patch_id = 
        patch_ids[n_hawk][random_patch];

    // error checking
#ifndef NDEBUG
    size_t n_hawk_check = 0;

    for (size_t breeder_i = 0; breeder_i < Npp; ++breeder_i)
    {
        n_hawk_check += MetaPop[random_patch_id].locals[breeder_i].hawk; 
    }

    assert(n_hawk_check == n_hawk);
#endif

    // fecundity selection absent
    // so whether we get immigrant 
    // only depends on dispersal probability
    size_t sampled = Npp + sample;

    // array to store in cumulative distribution over all sampled
    // individuals
    double choose_prob[sampled];

    // array to store patches from which the sampled individuals
    // come from
    size_t patch_origin[sampled];

    // array to store the ids of the individuals
    size_t individual_ids[sampled];

    // array with the individuals that may die
    size_t mortality_candidate_ids[Npp];
    size_t n_mortality_candidates = 0;

    double sumprob = 0;
    double d = 0;

    // at the same time calculate which local breeders
    // might die
    for (size_t local_i = 0; local_i < Npp; ++local_i)
    {
        // store the mortality candidates
        //
        // mortality hawk
        if (
                (
                 mortality_hawk && 
                    MetaPop[random_patch_id].locals[local_i].hawk 
                )
                ||
                (!mortality_hawk && 
                    !MetaPop[random_patch_id].locals[local_i].hawk 
                )

            )
        {
            mortality_candidate_ids[n_mortality_candidates++] = local_i;
        }
    }


    // now perform mortality and sample recruit

    // keep track of new state variables
    // to update patch type statistics and counters (see below)
    size_t n_hawk = n_hawk;

    // then mortality 
    assert(n_mortality_candidates > 0 && n_mortality_candidates <= Npp);

    // sample a random candidate for mortality
    size_t candidate_id = 
        mortality_candidate_ids[gsl_rng_uniform_int(r, n_mortality_candidates)];

    assert(candidate_id >= 0 && candidate_id < Npp);

    assert(
            (mortality_hawk &&
            MetaPop[random_patch_id].locals[candidate_id].hawk 
            )
            ||
            (!mortality_hawk &&
            !MetaPop[random_patch_id].locals[candidate_id].hawk 
            )
    );

    if (mortality_hawk)
    {
        --n_hawk_new;
    }


    // make new Kid
    Individual Kid;

    // sample local kid
    if (gsl_rng_uniform(r) < d)
    {
        create_kid(MetaPop[random_patch_id].locals[gsl_rng_uniform_int(r, Npp)],
                Individual);
    }
    else
    {
        // birth from a remote parent
        create_kid(MetaPop[gsl_rng_uniform_int(r, Npatches)].locals[
                gsl_rng_uniform_int(r, Npp)
                ], Individual);
    }

    // mortality and replacement
    MetaPop[random_patch_id].locals[candidate_id] = Kid;

    if (MetaPop[random_patch_id].locals[candidate_id].hawk)    
    {
        ++n_adapted_new;
    }

    // update statistics
    //
    // delete patch id in the corresponding stack
    // by replacing it with the last patch id in the stack
    // and then deleting the last element (by reducing the counter by 1)
    patch_ids[n_hawk][random_patch] =
        patch_ids[n_hawk][
            --Npatches[n_hawk] 
        ];

    // add patch id to the correct stack
    patch_ids[n_adapted_new][
        Npatches[n_adapted_new]++
    ] = random_patch_id;

    // we're done.
    // error checking
#ifndef NDEBUG
    n_hawk_check = 0;

    for (size_t breeder_i = 0; breeder_i < Npp; ++breeder_i)
    {
        n_hawk += MetaPop[random_patch_id].locals[breeder_i].hawk; 
    }

    assert(n_hawk_check == n_adapted_new);
#endif
}

// main function body
int main(int argc, char **argv)
{
    init_arguments(argc, argv);
    init_pop();
    write_parameters();
    write_data_headers();

    // there are this number of n probabilities:
    // 1. switching
    // 2. mortality adapted
    // 3. mortality maladapted adapted
    //
    // these can affect 2 * 3  different patch types
    // i. 2 x envt, ii. 3 x 0,1,2 adapted individuals per patch,
    // finally there are 
    size_t nprobs = 3 * 2 * 3;

    double probs[nprobs]; // vector with probabilities

    bool done; // variable to get out of nested loop
    double prob; // actual event probability sampled from 
                // the cumulative prob distribution

    // set up as a markovian process of numgen^2 timesteps
    // if we would have used a single variable then we would 
    // run into problems as the minimum specification of the maximum size
    // of size_t = ~65500
    for (generation = 0; generation < numgen; ++generation)
    {
        for (size_t generation2 = 0; generation2 < 80000; ++generation2)
        {
            // generate cumulative prob. distribution of possible
            // events
            double sum_probs = 0;
            size_t prob_count = 0;

            // loop through the different patch types 
            for (size_t patch_envt = 0; patch_envt < 2; ++patch_envt)
            {
                // number of adapted individuals in a patch
                for (size_t n_adapted = 0; n_adapted <= 2; ++n_adapted)
                {
                    // probability that a patch switches state
                    sum_probs += envt_switch[patch_envt] 
                        * Npatches[patch_envt][n_adapted];

                    probs[prob_count++] = sum_probs; 

                    // mortality adapted
                    sum_probs += (double) n_adapted / Npp 
                        * (n_adapted == 1 ? 
                             1.0 / (1.0 - C[patch_envt]) 
                             : 1.0 / (1.0 - c[patch_envt]) 
                          )
                        * Npatches[patch_envt][n_adapted];
                    
                    probs[prob_count++] = sum_probs; 

                    // mortality maladapted individual
                    sum_probs += (double) (Npp - n_adapted) / Npp 
                        * (n_adapted == 1 ? 
                             1.0 / (1.0 - C[patch_envt]) 
                             : 1.0 / (1.0 - c[patch_envt]) 
                          )
                        * Npatches[patch_envt][n_adapted];
                    
                    probs[prob_count++] = sum_probs; 
    
                    assert(prob_count <= nprobs);
                }
            }

            // decide which event will happen
            //
            // draw value from cumul distribution
            prob = gsl_rng_uniform(r) *  sum_probs;

            prob_count = 0;

            // variable to break out of all the loops
            done = false;

            // locate the event
            // in the cumulative distribution
            for (size_t patch_envt = 0; patch_envt < 2; ++patch_envt)
            {
                if (done)
                {
                    break;
                }

                for (size_t n_adapted = 0; n_adapted <= 2; ++n_adapted)
                {
                    if (prob <= probs[prob_count])
                    {
                        switch_patch_state(patch_envt, n_adapted);
                        done = true;
                        break;
                    }

                    ++prob_count;

                    // mortality adapted
                    if (prob <= probs[prob_count])
                    {
                        mortality(patch_envt, n_adapted, false);
                        done = true;
                        break;
                    }

                    ++prob_count;

                    // mortality maladapted
                    if (prob <= probs[prob_count])
                    {
                        mortality(patch_envt, n_adapted, true);
                        done = true;
                        break;
                    }
                    
                    ++prob_count;
                }
            } // end for patch_envt
    
            if (generation2 % skip2 == 0)
            {
                write_data();
            }
        } //end for size_t generation2
    } // end for size_t _generation

    write_data();

    exit(1);
}
