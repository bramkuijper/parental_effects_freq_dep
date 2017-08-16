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

size_t skip2 = numgen; // intervals of data output

// clutch size 
size_t clutch = 10;

// stats
size_t Nd = 0;
size_t Nh = 0;
size_t NdNext = 0;
size_t NhNext = 0;

// mortality of maladapted junior breeders in both envts
double c[2] = {0,0}; 

// mortality of adapted junior breeders in both envts
double C[2] = {0,0}; 

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
Individual NextGen[N*clutch];

// give the outputfile a unique name
string filename("sim_hd");
string filename_new(create_filename(filename));
ofstream DataFile(filename_new.c_str());  // output file 

// initialize the command line arguments to vary 
// the parameters
void init_arguments(int argc, char *argv[])
{
    c[0] = atof(argv[1]);
    c[1] = atof(argv[2]);
    C[0] = atof(argv[3]);
    C[1] = atof(argv[4]);
    mu_h = atof(argv[5]);
    sdmu_h = atof(argv[6]);
    init_pH = atof(argv[7]);
    init_pD = atof(argv[7]);
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

        Pop[i].is_hawk = gsl_rng_uniform(rng_r) < 0.5;

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

    assert(mother.pH[] >= 0);
    assert(mother.pH[1] <= 1.0);
    assert(mother.pD[1] >= 0);
    assert(mother.pD[1] <= 1.0);

    if (mother.is_hawk)
    {
        Kid.is_hawk = gsl_rng_uniform(r) < 
            0.5 * (mother.pH[0] + mother.pH[1]);
    }
    else
    {
        Kid.is_hawk = gsl_rng_uniform(r) < 
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

        // calculate payoffs
        if (Pop[ind1].is_hawk)
        {
            if (Pop[ind2].is_hawk) 
            {
                if (gsl_rng_uniform(rng_r) < 0.5)
                {
                    Pop[ind1].payoff = v;
                    Pop[ind2].payoff = c;
                }
                else
                {
                    Pop[ind1].payoff = c;
                    Pop[ind2].payoff = v;
                }
            }
            else
            {
                Pop[ind1].payoff = v;
                Pop[ind2].payoff = 0;
            }
        }
        else // ind1 is dove
        {
            if (Pop[ind2].is_hawk) 
            {
                Pop[ind1].payoff = 0;
                Pop[ind2].payoff = v;
            }
            else // both dove
            {
                Pop[ind1].payoff = .5 * v;
                Pop[ind2].payoff = .5 * v;
            }
        }

        // ok remove pairs
        --Npairstogo;

        // remove individual 1
        PopNext[Nnext++] = Pop[ind1];
        Pop[ind1] = Pop[--Ntot];

        // remove individual 2
        PopNext[Nnext++] = Pop[ind2];
        Pop[ind2] = Pop[--Ntot];

        assert(Nnext <= Ntot);
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
            if (PopNext[ind_j].cumul_payoff <= rand_val)
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
    double p1,p2,d,vard,varp1,varp2;
    double meanp1 = 0;
    double meanp2 = 0;
    double ssp1 = 0;
    double ssp2 = 0;

    double meand = 0;
    double ssd = 0;

    // loop through patches and individuals
    // and get stats on genotypes and patch freqs.
    for (size_t i = 0; i < Npatch; ++i)
    {
        for (size_t j = 0; j < Npp; ++j)
        {
            d = 0.5 * (MetaPop[i].locals[j].d[0] + MetaPop[i].locals[j].d[1]);
            p1 = 0.5 * (MetaPop[i].locals[j].p1[0] + MetaPop[i].locals[j].p1[1]);
            p2 = .5 * (MetaPop[i].locals[j].p2[0] + MetaPop[i].locals[j].p2[1]);

            meanp1+=p1;
            meanp2+=p2;
            meand+=d;
            ssp1+=p1*p1;
            ssp2+=p2*p2;
            ssd+=d*d;
        }    
    }

    DataFile << generation << ";"; 

    // write down all the patch frequency combinations
    for (size_t envt_i = 0; envt_i < 2; ++envt_i)
    {
        for (size_t n_adapted = 0; n_adapted <= 2; ++n_adapted)
        {
            DataFile << setprecision(5) << (double)Npatches[envt_i][n_adapted] / Npatch << ";";
        }
    }

    meanp1 /= Npatch * 2;
    meanp2 /= Npatch * 2;
    meand /= Npatch * 2;
    vard = ssd / (Npatch * 2) - meand*meand;
    varp1 = ssp1 / (Npatch * 2) - meanp1*meanp1;
    varp2 = ssp2 / (Npatch * 2) - meanp2*meanp2;

    DataFile 
        << setprecision(5) << meanp1 << ";"
        << setprecision(5) << meanp2 << ";"
        << setprecision(5) << varp1 << ";"
        << setprecision(5) << varp2 << ";"
        << setprecision(5) << meand << ";"
        << setprecision(5) << vard << endl;
}

// headers of the datafile
void write_data_headers()
{
    DataFile << "generation;";

    for (size_t envt_i = 0; envt_i < 2; ++envt_i)
    {
        for (size_t n_adapted = 0; n_adapted <= 2; ++n_adapted)
        {
            DataFile << setprecision(5) << "f_" << (envt_i+1) << "_" << (n_adapted + 1) << ";";
        }
    }
    
    DataFile << "meanp1;meanp2;varp1;varp2;meand;vard" << endl;
}

// parameters at the end of the sim
void write_parameters()
{
    total_time = time(NULL) - total_time;
    DataFile << endl << endl << "d0;" << d0 << endl
                << "type;" << "phenotype-dependent" << endl
                << "p0;" << h0 << endl
                << "mu_h;" << mu_h << endl
                << "sdmu_h;" << sdmu_h << endl
                << "numgen;" << numgen << endl
                << "k;" << k << endl
                << "s1;" << envt_switch[0] << endl
                << "s2;" << envt_switch[1] << endl
                << "c1;" << c[0] << endl
                << "c2;" << c[1] << endl
                << "C1;" << C[0] << endl
                << "C2;" << C[1] << endl
                << "seed;" << seed << endl
                << "Npatches;" << Npatches << endl
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
    } 

    write_data();

    exit(1);
}
