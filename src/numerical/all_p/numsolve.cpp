#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <cassert>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include "bramauxiliary.h"


using namespace std;

string filename("iter_freqdep");
string filename_new(create_filename(filename));
ofstream DataFile(filename_new.c_str());  // output file 

struct rparams
{
    // VARS
};

double bound0(double val)
{
    if (val < 0)
    {
        val = 0.001;
    }

    return(val);
}

double bound01(double val)
{
    val = val < 0 ? 0.0001 : val > 1.0 ? 0.9999 : val;

    return(val);
}

// recursions of all the patch frequencies
int psys_recur(void *params, gsl_vector *f)
{
    //VARFUNC
    
    
    // PATCHRECUR
    
    return GSL_SUCCESS;
}

void relcoeff(void *params, gsl_vector *f)
{
    //VARFUNC


    // RELATEDNESS
}


void reproductive_values(void *params, gsl_vector *f)
{
    // VARFUNC


    // REPVALS
}

void selgrads(void *params, gsl_vector *f)
{
    // VARFUNC

    // SELGRADS

}


void write_params(void *params)
{
    // VARFUNC


    // WRITEPARAMS
}


void write_data(void *params, int time)
{
    // VARFUNC


    if (time < 0)
    {
        // HEADERWRT
    }

    // WRITEDATA
}


int main (int argc, char **argv)
{
    int max_iter = atoi(argv[1]);

    // initialize the vectors that contain the variables
    // functions solve for
    //
    // vector for patch frequencies
    gsl_vector *x_p = gsl_vector_alloc(2);
    
    // vector for relatedness
    gsl_vector *x_r = gsl_vector_alloc(3);

    // vector for reproductive values
    gsl_vector *x_v = gsl_vector_alloc(4);

    // vector for selection gradients
    gsl_vector *x_selgrad = gsl_vector_alloc(4);

    // initialize the struct with parameters
    struct rparams paramstruct; 

    // initialize command line argument things
    // see generate_cpp.py 
    // ARGINIT
   

    
    
    // ranges for cycling: sometimes numerical iterations
    // won't resolve as solutions slightly cycle around
    // the final value. To see whether cyling behaviour
    // is repetitive over time, however, we create these
    // vectors that stores a series of past values of 
    // the values for the switching rates and compares
    // those to values that are currently found
    gsl_vector *phh1_range = gsl_vector_alloc(10);
    gsl_vector *phh2_range = gsl_vector_alloc(10);
    gsl_vector *pdh1_range = gsl_vector_alloc(10);
    gsl_vector *pdh0_range = gsl_vector_alloc(10);

    for (int ik = 0; ik < 10; ++ik)
    {
        // initialize the vector
        gsl_vector_set(phh1_range, ik, 0);        
        gsl_vector_set(phh2_range, ik, 0);        
        gsl_vector_set(pdh1_range, ik, 0);        
        gsl_vector_set(pdh0_range, ik, 0);        
    }

    // write the initial data set
    write_data(&paramstruct,-1);

    // iterate
    int iter;
    for (iter = 0; iter < max_iter ; ++iter)
    {
        // patch frequencies
        psys_recur(&paramstruct, x_p);

        paramstruct.f_0= gsl_vector_get(x_p,0);
        paramstruct.f_1= gsl_vector_get(x_p,1);

        assert(isnan(paramstruct.f_0) == 0);
        assert(isnan(paramstruct.f_1) == 0);

        // relatedness
        relcoeff(&paramstruct, x_r);        
        
        paramstruct.rhh = gsl_vector_get(x_r, 0);
        paramstruct.rhd = gsl_vector_get(x_r, 1);
        paramstruct.rdd = gsl_vector_get(x_r, 2);
        
        assert(isnan(paramstruct.rhh) == 0);
        assert(isnan(paramstruct.rhd) == 0);
        assert(isnan(paramstruct.rdd) == 0);

        // reproductive values
        reproductive_values(&paramstruct, x_v);
        
        paramstruct.vh_1 = gsl_vector_get(x_v, 0);
        paramstruct.vh_2 = gsl_vector_get(x_v, 1);
        paramstruct.vd_0 = gsl_vector_get(x_v, 2);
        paramstruct.vd_1 = gsl_vector_get(x_v, 3);
        
        assert(isnan(paramstruct.vh_1) == 0);
        assert(isnan(paramstruct.vh_2) == 0);
        assert(isnan(paramstruct.vd_0) == 0);
        assert(isnan(paramstruct.vd_1) == 0);

        // selection gradients
        selgrads(&paramstruct, x_selgrad);

        bool condition_phh1 = (fabs(paramstruct.phh1 - gsl_vector_get(x_selgrad, 0)) < 1e-7 || paramstruct.phh1 >= 0.999) || paramstruct.phh1 <= 0.001;
        bool condition_phh2 = (fabs(paramstruct.phh2 - gsl_vector_get(x_selgrad, 1)) < 1e-7 || paramstruct.phh2 >= 0.999) || paramstruct.phh2 <= 0.001;
        bool condition_pdh1 = (fabs(paramstruct.pdh1 - gsl_vector_get(x_selgrad, 2)) < 1e-7 || paramstruct.pdh1 >= 0.999) || paramstruct.pdh1 <= 0.001;
        bool condition_pdh0 = (fabs(paramstruct.pdh0 - gsl_vector_get(x_selgrad, 3)) < 1e-7 || paramstruct.pdh0 >= 0.999) || paramstruct.pdh0 <= 0.001;

        if (condition_phh2 && 
                condition_phh1 && 
                condition_pdh1 && 
                condition_pdh0)
        {
            paramstruct.phh1 = gsl_vector_get(x_selgrad, 0);
            paramstruct.phh2 = gsl_vector_get(x_selgrad, 1);
            paramstruct.pdh1 = gsl_vector_get(x_selgrad, 2);
            paramstruct.pdh0 = gsl_vector_get(x_selgrad, 3);

            break;
        }
        paramstruct.phh1 = gsl_vector_get(x_selgrad, 0);
        paramstruct.phh2 = gsl_vector_get(x_selgrad, 1);
        paramstruct.pdh1 = gsl_vector_get(x_selgrad, 2);
        paramstruct.pdh0 = gsl_vector_get(x_selgrad, 3);
        
        assert(isnan(paramstruct.phh1) == 0);
        assert(isnan(paramstruct.phh2) == 0);
        assert(isnan(paramstruct.pdh1) == 0);
        assert(isnan(paramstruct.pdh0) == 0);
       
        // check for stable nonequilibrium behaviour 
        if (iter > 50000)
        {
            bool found_in_range1 = false;
            bool found_in_range2 = false;
            bool found_in_range3 = false;
            bool found_in_range4 = false;

            bool done = false;

            for (int ik = 0; ik < 10; ++ik)
            {
                if (!found_in_range1 && fabs(gsl_vector_get(phh1_range, ik)-paramstruct.phh1) < 1e-10)
                {
                    found_in_range1 = true;
                }
                if (!found_in_range2 && fabs(gsl_vector_get(phh2_range, ik)-paramstruct.phh2) < 1e-10)
                {
                    found_in_range2 = true;
                }
                if (!found_in_range3 && fabs(gsl_vector_get(pdh1_range, ik)-paramstruct.pdh1) < 1e-10)
                {
                    found_in_range3 = true;
                }
                if (!found_in_range4 && fabs(gsl_vector_get(pdh0_range, ik)-paramstruct.pdh0) < 1e-10)
                {
                    found_in_range4 = true;
                }

                if (found_in_range1 && found_in_range2 && found_in_range3 && found_in_range4)
                {
                    done = true;
                    break;
                }
            }

            if (done)
            {
                break;
            }
        }


        for (int ik = 9; ik > 0; --ik)
        {
            gsl_vector_set(phh1_range, ik, gsl_vector_get(phh1_range, ik - 1));
            gsl_vector_set(phh2_range, ik, gsl_vector_get(phh2_range, ik - 1));
            gsl_vector_set(pdh1_range, ik, gsl_vector_get(pdh1_range, ik - 1));
            gsl_vector_set(pdh0_range, ik, gsl_vector_get(pdh0_range, ik - 1));
        }

        gsl_vector_set(phh1_range, 0, paramstruct.phh1);
        gsl_vector_set(phh2_range, 0, paramstruct.phh2);
        gsl_vector_set(pdh1_range, 0, paramstruct.pdh1);
        gsl_vector_set(pdh0_range, 0, paramstruct.pdh0);

        if (iter % 100 == 0)
        {
            write_data(&paramstruct,iter);
        }
    }

    write_data(&paramstruct,iter);
    write_params(&paramstruct);


    gsl_vector_free(x_p);
    gsl_vector_free(x_r);
    gsl_vector_free(x_v);
    gsl_vector_free(x_selgrad);
}
