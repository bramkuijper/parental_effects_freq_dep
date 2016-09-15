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
unsigned long int NumGen = 1e10; // maximum number of generations
int sample = 20;

double h0 = 0.5;

int seed = -1; // the random seed
time_t total_time;  // track the total time the simulation is running
unsigned long int generation = 0; // track the current generation 

// output once every 10^6 generations
unsigned long int skip = 1e6;

double d = 0.1; // dispersal probability
double v = 0.1;
double c = 0.1;


// mutation rates
double mu_h = 0.01;
double sdmu_h = 0.01;

struct Individual {

    double phh1[2];
    double phh2[2];
    double pdh0[2];
    double pdh1[2];

    bool hawk;

};

struct Patch {
    Individual locals[Npp];
    size_t nhawk;
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
    h0 = atof(argv[4]);
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
                ++n_hawk;
            }
            else
            {
                // randomly assign initial phenotypes
                MetaPop[i].locals[j].hawk = false;
            }

            // initialize allelic values
            for (int allele = 0; allele < 2; ++allele)
            {
                MetaPop[i].locals[j].phh1[allele] = h0;
                MetaPop[i].locals[j].phh2[allele] = h0;
                MetaPop[i].locals[j].pdh0[allele] = h0;
                MetaPop[i].locals[j].pdh1[allele] = h0;
            }
        }

        // increase the counter of patches 
        // containing n_hawk individuals
        //
        size_t Npatches_of_this_type = 
            Npatches[n_hawk];

        MetaPop[i].nhawk = n_hawk;

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
void create_kid(Individual &mother, Individual &Kid, size_t nhawk)
{
    double expr_h = -1.0;

    switch(nhawk)
    {
        case 0:
            expr_h = 0.5*(mother.pdh0[0] + mother.pdh0[1]);
            break;
        case 1:
            expr_h = mother.hawk ? 
                (0.5 * (mother.phh1[0] + mother.phh1[1]))
                :
                (0.5 * (mother.pdh1[0] + mother.pdh1[1]));
            break;
        case 2:
            expr_h = 0.5 * (mother.phh2[0] + mother.phh2[1]);
            break;
        default:
            assert(nhawk <= 2);
            break;
    }

    assert(expr_h >= 0 && expr_h <= 1.0);
    
    // mother's phenotype is one

    // offspring gets hawk or dove phenotype
    Kid.hawk = gsl_rng_uniform(r) < expr_h;

    // inherit alleles
    for (int allele_i = 0; allele_i < 2; ++allele_i)
    {
        Kid.phh1[allele_i] = mut_h(mother.phh1[allele_i]);
        Kid.phh2[allele_i] = mut_h(mother.phh2[allele_i]);

        Kid.pdh0[allele_i] = mut_h(mother.pdh0[allele_i]);
        Kid.pdh1[allele_i] = mut_h(mother.pdh1[allele_i]);
    }
}

// do the stats on the data
void write_data()
{
    // get variance and means
    double phh2, phh1, pdh0, pdh1;
    double varphh2, varphh1, varpdh0, varpdh1;
    double meanphh2, meanphh1, meanpdh0, meanpdh1;
    double ssphh2, ssphh1, sspdh0, sspdh1;

    meanphh2 = 0;
    meanphh1 = 0;
    meanpdh0 = 0;
    meanpdh1 = 0;

    // loop through patches and individuals
    // and get stats on genotypes and patch freqs.
    for (size_t i = 0; i < Npatch; ++i)
    {
        for (size_t j = 0; j < Npp; ++j)
        {
            phh1 = 0.5 * (MetaPop[i].locals[j].phh1[0] + MetaPop[i].locals[j].phh1[1]);
            phh2 = 0.5 * (MetaPop[i].locals[j].phh2[0] + MetaPop[i].locals[j].phh2[1]);
            pdh0 = .5 * (MetaPop[i].locals[j].pdh0[0] + MetaPop[i].locals[j].pdh0[1]);
            pdh1 = .5 * (MetaPop[i].locals[j].pdh1[0] + MetaPop[i].locals[j].pdh1[1]);

            meanphh1+=phh1;
            meanphh2+=phh2;
            meanpdh0+=pdh0;
            meanpdh1+=pdh1;

            ssphh1+=phh1*phh1;
            ssphh2+=phh2*phh2;
            sspdh0+=pdh0*pdh0;
            sspdh1+=pdh1*pdh1;
        }    
    }

    DataFile << generation << ";"; 

    // write down all the patch frequency combinations
    for (size_t n_hawk = 0; n_hawk <= 2; ++n_hawk)
    {
        DataFile << setprecision(5) 
            << (double)Npatches[n_hawk] / Npatch << ";";
    }

    meanphh1 /= Npatch * 2;
    meanphh2 /= Npatch * 2;
    meanpdh0 /= Npatch * 2;
    meanpdh1 /= Npatch * 2;
    varphh2 = ssphh2 / (Npatch * 2) - meanphh2*meanphh2;
    varphh1 = ssphh1 / (Npatch * 2) - meanphh1*meanphh1;
    varpdh0 = sspdh0 / (Npatch * 2) - meanpdh0*meanpdh0;
    varpdh1 = sspdh1 / (Npatch * 2) - meanpdh1*meanpdh1;

    DataFile 
        << setprecision(5) << meanphh1 << ";"
        << setprecision(5) << meanphh2 << ";"
        << setprecision(5) << meanpdh1 << ";"
        << setprecision(5) << meanpdh0 << ";"
        << setprecision(5) << varphh2 << ";"
        << setprecision(5) << varphh1 << ";"
        << setprecision(5) << varpdh0 << ";" 
        << setprecision(5) << varpdh1 << ";" << endl;
}

// headers of the datafile
void write_data_headers()
{
    DataFile << "generation;";

    for (size_t n_hawk = 0; n_hawk <= 2; ++n_hawk)
    {
        DataFile << setprecision(5) << "f_" << (n_hawk) << ";";
    }
    
    DataFile << "meanphh1;meanphh2;meanpdh1;meanpdh0;varphh2;varphh1;varpdh0;varpdh1;" << endl;
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
                << "Npatches;" << Npatch << endl
                << "runtime;" << total_time << endl;
}


// mortality of an individual in one of the patches 
//
// size_t n_hawk: a patch containing n_hawk adapted individuals
// bool mortality_maladapted: true, maladapted guy dies; false, adapted guy dies;

void mortality(size_t n_hawk, bool mortality_hawk)
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

    // array with the individuals that may die
    size_t mortality_candidate_ids[Npp];
    size_t n_mortality_candidates = 0;

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
    size_t n_hawk_new = n_hawk;

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

    size_t random_remote_patch_id = 0;

    // sample local kid
    if (gsl_rng_uniform(r) < d)
    {
        create_kid(
                MetaPop[random_patch_id].locals[gsl_rng_uniform_int(r, Npp)],
                Kid,
                n_hawk
                );
    }
    else
    {
        // birth from a remote parent
        random_remote_patch_id = gsl_rng_uniform_int(r, Npatch);

        create_kid(
                MetaPop[random_remote_patch_id].locals[gsl_rng_uniform_int(r, Npp)], 
                Kid,
                MetaPop[random_remote_patch_id].nhawk
                );
    }

    // mortality and replacement
    MetaPop[random_patch_id].locals[candidate_id] = Kid;

    // if newborn is hawk, increment count
    if (MetaPop[random_patch_id].locals[candidate_id].hawk)    
    {
        ++n_hawk_new;
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
    patch_ids[n_hawk_new][
        Npatches[n_hawk_new]++
    ] = random_patch_id;

    MetaPop[random_patch_id].nhawk = n_hawk_new;

    // we're done.
    // error checking
#ifndef NDEBUG
    n_hawk_check = 0;

    for (size_t breeder_i = 0; breeder_i < Npp; ++breeder_i)
    {
        n_hawk_check += MetaPop[random_patch_id].locals[breeder_i].hawk; 
    }

    assert(n_hawk_check == n_hawk_new);
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
    // 1. mortality hawk on a 1 hawk 1 dove patch
    // 2. mortality hawk on a 2 hawk patch
    // 3. mortality dove on a 1 hawk 1 dove patch
    // 4. mortality dove on a 2 dove patch
    //
    size_t nprobs = 4;

    double probs[nprobs]; // vector with probabilities

    double prob; // actual event probability sampled from 
                // the cumulative prob distribution

    // run the markov process
    for (generation = 0; generation < NumGen; ++generation)
    {
        // generate cumulative prob. distribution of possible
        // events
        double sum_probs = 0;

        /// 1. mortality hawk on a 1 hawk 1 dove patch
        sum_probs += .5 * (1.0-v) * Npatches[1];
        probs[0] = sum_probs;

        // 2. mortality hawk on a 2 hawk patch
        sum_probs += (1.0-(.5*v - .5*c)) * Npatches[2];
        probs[1] = sum_probs;

        // 3. mortality dove on a 1 hawk 1 dove patch
        sum_probs += .5 * (1.0-0) * Npatches[1];
        probs[2] = sum_probs;
    
        // 4. mortality dove on a 2 dove patch
        sum_probs += (1.0-.5*v) * Npatches[0];
        probs[3] = sum_probs;

        // decide which event will happen
        //
        // draw value from cumul distribution
        prob = gsl_rng_uniform(r) *  sum_probs;

        // 1. mortality hawk on a 1 hawk 1 dove patch
        //
        if (prob <= probs[0])
        {
            mortality(1, true);
        } else if (prob <= probs[1])
        {
            mortality(2, true);
        }
        else if (prob <= probs[2])
        {
            mortality(1, false);
        }
        else if (prob <= probs[3])
        {
            mortality(0, false);
        }
    
        if (generation % skip == 0)
        {
            write_data();
        }
    } // end for size_t _generation

    write_data();

    exit(1);
}
