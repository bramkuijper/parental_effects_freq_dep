// Discrete-time hawk-dove game
// with the evolution of parental effects

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
gsl_rng *rng_r; // gnu scientific rng 

const size_t N = 5000; 
const size_t numgen = 40000;

int seed = -1; // the random seed
time_t total_time;  // track the total time the simulation is running


size_t generation = 0; // track the current generation (for debugging only)

size_t skip = 10; // intervals of data output

// clutch size 
size_t clutch = 10;

// stats
size_t Nd = 0;
size_t Nh = 0;
size_t NdNext = 0;
size_t NhNext = 0;

// v, c
double v = 0;
double c = 0;

// mutation rates
double mu_h = 0.01;
double sdmu_h = 0.01;

double init_pH = 0.5;
double init_pD = 0.5;

// diploid individuals
struct Individual
{
    double pH[2]; // proportion z1 offspring by z1 junior
    double pD[2]; // proportion z1 offspring by z2 junior

    bool is_hawk; 

    double payoff; // absolute payoff
    double cumul_payoff; // position in cumulative distribution of payoffs
};

// definie populations
Individual Pop[N];
Individual PopNext[N];

// give the outputfile a unique name
string filename("sim_hd");
string filename_new(create_filename(filename));
ofstream DataFile(filename_new.c_str());  // output file 

// initialize the command line arguments to vary 
// the parameters
void init_arguments(int argc, char *argv[])
{
    v = atof(argv[1]);
    c = atof(argv[2]);
    mu_h = atof(argv[3]);
    sdmu_h = atof(argv[4]);
    init_pH = atof(argv[5]);
    init_pD = atof(argv[6]);
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
    rng_r = gsl_rng_alloc(T);
    gsl_rng_set(rng_r, seed);

    Nh = 0;
    Nd = 0;

    // initialize all the patches
    for (size_t i = 0; i < N; ++i)
    {
        for (size_t allele_i = 0; 
                allele_i < 2;
                ++allele_i)
        {
            Pop[i].pH[allele_i] = init_pH;
            Pop[i].pD[allele_i] = init_pD;
        }

        Pop[i].is_hawk = gsl_rng_uniform(rng_r) < v/c;

        Pop[i].payoff = 0;
        Pop[i].cumul_payoff = 0;
    }
}

double mut_h(double p)
{
    assert(p >= 0 && p <= 1);
    p += gsl_rng_uniform(rng_r) < mu_h ? gsl_ran_gaussian(rng_r,sdmu_h) : 0;
    p = p < 0 ? 0 : p > 1 ? 1 : p;

    return(p);
}

// make a kid out of parent asexually
// add it either to local offspring pool
// or to pool of dispersers
void create_kid(Individual &mother, Individual &Kid)
{

    assert(mother.pH[0] >= 0);
    assert(mother.pH[0] <= 1.0);
    assert(mother.pD[0] >= 0);
    assert(mother.pD[0] <= 1.0);

    assert(mother.pH[1] >= 0);
    assert(mother.pH[1] <= 1.0);
    assert(mother.pD[1] >= 0);
    assert(mother.pD[1] <= 1.0);

    if (mother.is_hawk)
    {
        Kid.is_hawk = gsl_rng_uniform(rng_r) < 
            0.5 * (mother.pH[0] + mother.pH[1]);
    }
    else
    {
        Kid.is_hawk = gsl_rng_uniform(rng_r) < 
            1.0 - 0.5 * (mother.pD[0] + mother.pD[1]);
    }

    // inherit alleles
    for (int allele_i = 0; allele_i < 2; ++allele_i)
    {
        Kid.pH[allele_i] = mut_h(mother.pH[allele_i]);
        Kid.pD[allele_i] = mut_h(mother.pD[allele_i]);
    }

    Kid.payoff = 0;
    Kid.cumul_payoff = 0;
}

// let hawks, doves interact among one another
// then do survival
void survive()
{
    size_t Nnext = 0;

    size_t ind1, ind2;

    size_t Npairstogo = N/2;

    size_t Ntot = N;

    // Let N/2 
    for (size_t pair_i = 0; pair_i < Npairstogo; ++pair_i) 
    {
        // no pairs left
        if (Ntot < 2)
        {
            break;
        }

        ind1 = gsl_rng_uniform_int(rng_r, Ntot);

        do {
            ind2 = gsl_rng_uniform_int(rng_r, Ntot);
        }
        while (ind1 == ind2);

        Pop[ind1].payoff = c;
        Pop[ind2].payoff = c;

        // calculate payoffs
        if (Pop[ind1].is_hawk)
        {
            if (Pop[ind2].is_hawk) 
            {
                if (gsl_rng_uniform(rng_r) < 0.5)
                {
                    Pop[ind1].payoff += v;
                    Pop[ind2].payoff += -c;
                }
                else
                {
                    Pop[ind1].payoff += -c;
                    Pop[ind2].payoff += v;
                }
            }
            else
            {
                Pop[ind1].payoff += v;
                Pop[ind2].payoff += 0;
            }
        }
        else // ind1 is dove
        {
            if (Pop[ind2].is_hawk) 
            {
                Pop[ind1].payoff += 0;
                Pop[ind2].payoff += v;
            }
            else // both dove
            {
                Pop[ind1].payoff += .5 * v;
                Pop[ind2].payoff += .5 * v;
            }
        }

        // remove individual 1
        PopNext[Nnext++] = Pop[ind1];
        Pop[ind1] = Pop[--Ntot];

        // remove individual 2
        PopNext[Nnext++] = Pop[ind2];
        Pop[ind2] = Pop[--Ntot];

//        cout << Nnext << " " << pair_i << Ntot << endl;

        assert(Nnext <= N);
    }

    double sum_payoffs = 0;

    // make cumulative distribution of payoffs
    for (size_t ind_i = 0; ind_i < Nnext; ++ind_i)
    {
        sum_payoffs = PopNext[ind_i].cumul_payoff = 
            sum_payoffs + PopNext[ind_i].payoff;
    }

    // now make the next generation
    for (size_t ind_i = 0; ind_i < N; ++ind_i)
    {
        // select random number
        double rand_val = gsl_rng_uniform(rng_r) 
            * sum_payoffs;

        // select from cumulative distribution
        for (size_t ind_j = 0; ind_j < N; ++ind_j)
        {
            if (rand_val <= PopNext[ind_j].cumul_payoff)
            {
                Individual Kid;
                create_kid(PopNext[ind_j], Kid);

                Pop[ind_i] = Kid;

                break;
            }
        }
    }
}


// do the stats on the data
void write_data()
{
    // get variance and means
    double pH, pD;
    double meanpH = 0;
    double meanpD = 0;
    double sspH = 0;
    double sspD = 0;
    double freq_hawk = 0;
    double varpH, varpD;

    // loop through all individuals
    // and get stats on genotypes and patch freqs.
    for (size_t ind_i = 0; ind_i < N; ++ind_i)
    {
        pH = 0.5 * (Pop[ind_i].pH[0] + Pop[ind_i].pH[1]);
        pD = 0.5 * (Pop[ind_i].pD[0] + Pop[ind_i].pD[1]);

        freq_hawk += Pop[ind_i].is_hawk;

        meanpH+=pH;
        meanpD+=pD;
        sspH+=pH*pH;
        sspD+=pD*pD;
    }

    DataFile << generation << ";"; 

    meanpH /= N;
    meanpD /= N;
    varpH = sspH / N - meanpH * meanpH;
    varpD = sspD / N - meanpD * meanpD;

    freq_hawk /= N;

    DataFile 
        << setprecision(5) << meanpH << ";"
        << setprecision(5) << meanpD << ";"
        << setprecision(5) << varpH << ";"
        << setprecision(5) << varpD << ";" 
        << setprecision(5) << freq_hawk << ";" 
        << endl;
}

// headers of the datafile
void write_data_headers()
{
    DataFile << "generation;meanpH;meanpD;varpH;varpD;freq_hawk;" << endl;
}

// parameters at the end of the sim
void write_parameters()
{
    total_time = time(NULL) - total_time;
    DataFile << endl << endl 
                << "v;" << v << endl
                << "c;" << c << endl
                << "init_pH;" << init_pH << endl
                << "init_pD;" << init_pD << endl
                << "mu_h;" << mu_h << endl
                << "sdmu_h;" << sdmu_h << endl
                << "numgen;" << numgen << endl
                << "seed;" << seed << endl
                << "N;" << N << endl
                << "runtime;" << total_time << endl;
}


// main function body
int main(int argc, char **argv)
{
    init_arguments(argc, argv);
    init_pop();
    write_parameters();
    write_data_headers();

    for (generation = 0; generation < numgen; ++generation)
    {
        survive();

        if (generation % skip == 0)
        {
            write_data();
        }
    } 

    write_data();

    exit(1);
}
