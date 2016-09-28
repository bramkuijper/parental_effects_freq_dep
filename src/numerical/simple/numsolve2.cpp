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
    // 
double d;
double md_0;
double mh_1;
double md_1;
double mh_2;
double pdhinit;
double phhinit;
double f_0;
double f_1;
double pdh;
double phh;
double rhh;
double rhd;
double rdd;
double vh_1;
double vh_2;
double vd_0;
double vd_1;

};

double bound(double val)
{
    val = val < 0 ? 0.0001 : val > 1.0 ? 0.9999 : val;

    return(val);
}


// recursions of all the patch frequencies
int psys_recur(void *params, gsl_vector *f)
{
    //
double d = ((struct rparams *) params)->d;
double md_0 = ((struct rparams *) params)->md_0;
double mh_1 = ((struct rparams *) params)->mh_1;
double md_1 = ((struct rparams *) params)->md_1;
double mh_2 = ((struct rparams *) params)->mh_2;
double pdhinit = ((struct rparams *) params)->pdhinit;
double phhinit = ((struct rparams *) params)->phhinit;
double f_0 = ((struct rparams *) params)->f_0;
double f_1 = ((struct rparams *) params)->f_1;
double pdh = ((struct rparams *) params)->pdh;
double phh = ((struct rparams *) params)->phh;
double rhh = ((struct rparams *) params)->rhh;
double rhd = ((struct rparams *) params)->rhd;
double rdd = ((struct rparams *) params)->rdd;
double vh_1 = ((struct rparams *) params)->vh_1;
double vh_2 = ((struct rparams *) params)->vh_2;
double vd_0 = ((struct rparams *) params)->vd_0;
double vd_1 = ((struct rparams *) params)->vd_1;

    
    
    // 

double df0dt = 0;

double df1dt = 0;

for (int rtime=0;rtime<1e7;++rtime)
{

 df0dt=f_0 + 0.01*(-(f_0*(2*(1 - d)*pdh + d*(2*pdh*f_0 + 2*phh*(1 - f_0 - f_1) + (pdh + phh)*f_1))*md_0) + f_1*(1 + (-((1 - d)*(pdh + phh)) - d*(2*pdh*f_0 + 2*phh*(1 - f_0 - f_1) + (pdh + phh)*f_1))/2.)*mh_1);

df1dt=f_1 + 0.01*(f_0*(2*(1 - d)*pdh + d*(2*pdh*f_0 + 2*phh*(1 - f_0 - f_1) + (pdh + phh)*f_1))*md_0 - (f_1*((1 - d)*(pdh + phh) + d*(2*pdh*f_0 + 2*phh*(1 - f_0 - f_1) + (pdh + phh)*f_1))*md_1)/2. - f_1*(1 + (-((1 - d)*(pdh + phh)) - d*(2*pdh*f_0 + 2*phh*(1 - f_0 - f_1) + (pdh + phh)*f_1))/2.)*mh_1 + 2*(1 - f_0 - f_1)*(1 + (-2*(1 - d)*phh - d*(2*pdh*f_0 + 2*phh*(1 - f_0 - f_1) + (pdh + phh)*f_1))/2.)*mh_2);

if (fabs(df0dt-f_0)<1e-7 && fabs(df1dt-f_1)<1e-7)
{

f_0=bound(df0dt);

f_1=bound(df1dt);

break;
}
f_0=bound(df0dt);

f_1=bound(df1dt);

}
gsl_vector_set(f,0,df0dt);
gsl_vector_set(f,1,df1dt);
    
    return GSL_SUCCESS;
}

void relcoeff(void *params, gsl_vector *f)
{
    //
double d = ((struct rparams *) params)->d;
double md_0 = ((struct rparams *) params)->md_0;
double mh_1 = ((struct rparams *) params)->mh_1;
double md_1 = ((struct rparams *) params)->md_1;
double mh_2 = ((struct rparams *) params)->mh_2;
double pdhinit = ((struct rparams *) params)->pdhinit;
double phhinit = ((struct rparams *) params)->phhinit;
double f_0 = ((struct rparams *) params)->f_0;
double f_1 = ((struct rparams *) params)->f_1;
double pdh = ((struct rparams *) params)->pdh;
double phh = ((struct rparams *) params)->phh;
double rhh = ((struct rparams *) params)->rhh;
double rhd = ((struct rparams *) params)->rhd;
double rdd = ((struct rparams *) params)->rdd;
double vh_1 = ((struct rparams *) params)->vh_1;
double vh_2 = ((struct rparams *) params)->vh_2;
double vd_0 = ((struct rparams *) params)->vd_0;
double vd_1 = ((struct rparams *) params)->vd_1;



    // 
double rhhtplus1 = -(((-1 + d)*(-((-1 + d)*(pdh - phh)*pow(f_1,2)*md_1*((2 - (-1 + d)*(-1 + pdh)*f_0*md_0)*md_1 + (2*(pdh + phh) + (-1 + d)*phh*f_0*md_0)*mh_1)) + 2*phh*(1 - f_0 - f_1)*(2 - (-1 + d)*(-1 + pdh)*f_0*md_0)*mh_2 + f_1*(2*(2*phh - (-1 + d)*(-phh + pdh*(pdh + phh))*f_0*md_0)*md_1 + (-1 + d)*(1 - f_0 - f_1)*((pdh - phh)*(-2 + (-1 + d)*(-1 + pdh)*f_0*md_0)*md_1 + phh*(2*phh - (-1 + d)*(pdh - phh)*f_0*md_0)*mh_1)*mh_2)))/(2*(4 - 2*phh*(1 - f_0 - f_1)*mh_2 + 2*d*phh*(1 - f_0 - f_1)*mh_2 - (-1 + d)*f_1*(phh*mh_1*(-2 - (-1 + d)*phh*(1 - f_0 - f_1)*mh_2) + md_1*(-2 - (-1 + d)*phh*(1 - f_0 - f_1)*mh_2 + pdh*(2 + (-1 + d)*(1 - f_0 - f_1)*mh_2)))) + (-1 + d)*f_0*md_0*(-2*(-1 + pdh)*(2 + (-1 + d)*phh*(1 - f_0 - f_1)*mh_2) + (-1 + d)*f_1*(-((pdh - phh)*mh_1*(2 + (-1 + d)*phh*(1 - f_0 - f_1)*mh_2)) + (-1 + pdh)*md_1*(-2 - (-1 + d)*phh*(1 - f_0 - f_1)*mh_2 + pdh*(2 + (-1 + d)*(1 - f_0 - f_1)*mh_2))))));

double rhdtplus1 = -(((-1 + d)*(2*(2*pdh*f_0*md_0 + (1 - f_0 - f_1)*(2 - 2*phh + (-1 + d)*(1 - phh + pdh*(-1 + 2*phh))*f_0*md_0)*mh_2) + f_1*(2*pdh*mh_1*(2 + (-1 + d)*phh*(1 - f_0 - f_1)*mh_2) + (-2 + (-1 + d)*(-1 + pdh)*f_0*md_0)*md_1*(-2 - (-1 + d)*phh*(1 - f_0 - f_1)*mh_2 + phh*(2 + (-1 + d)*(1 - f_0 - f_1)*mh_2)))))/(2*(4 - 2*phh*(1 - f_0 - f_1)*mh_2 + 2*d*phh*(1 - f_0 - f_1)*mh_2 - (-1 + d)*f_1*(phh*mh_1*(-2 - (-1 + d)*phh*(1 - f_0 - f_1)*mh_2) + md_1*(-2 - (-1 + d)*phh*(1 - f_0 - f_1)*mh_2 + pdh*(2 + (-1 + d)*(1 - f_0 - f_1)*mh_2)))) + (-1 + d)*f_0*md_0*(-2*(-1 + pdh)*(2 + (-1 + d)*phh*(1 - f_0 - f_1)*mh_2) + (-1 + d)*f_1*(-((pdh - phh)*mh_1*(2 + (-1 + d)*phh*(1 - f_0 - f_1)*mh_2)) + (-1 + pdh)*md_1*(-2 - (-1 + d)*phh*(1 - f_0 - f_1)*mh_2 + pdh*(2 + (-1 + d)*(1 - f_0 - f_1)*mh_2))))));

double rddtplus1 = -(((-1 + d)*(f_1*mh_1*(2*(2 + (-1 + d)*pow(pdh,2)*f_1*md_1 + (-1 + d)*phh*f_1*(-((-2 + phh)*md_1) + mh_1) + pdh*(-2 - (-1 + d)*f_1*(2*md_1 + mh_1))) + (-1 + d)*(1 - f_0 - f_1)*(-2*(1 + (-2 + pdh)*phh + (-1 + phh)*phh) + (-1 + d)*(pdh - phh)*f_1*((-1 + pdh)*md_1 - phh*mh_1))*mh_2) + f_0*md_0*(-2*(-1 + pdh)*(2 + (-1 + d)*phh*(1 - f_0 - f_1)*mh_2) + (-1 + d)*f_1*(-((pdh - phh)*mh_1*(2 + (-1 + d)*phh*(1 - f_0 - f_1)*mh_2)) + (-1 + pdh)*md_1*(-2 - (-1 + d)*phh*(1 - f_0 - f_1)*mh_2 + pdh*(2 + (-1 + d)*(1 - f_0 - f_1)*mh_2))))))/(2*(4 - 2*phh*(1 - f_0 - f_1)*mh_2 + 2*d*phh*(1 - f_0 - f_1)*mh_2 - (-1 + d)*f_1*(phh*mh_1*(-2 - (-1 + d)*phh*(1 - f_0 - f_1)*mh_2) + md_1*(-2 - (-1 + d)*phh*(1 - f_0 - f_1)*mh_2 + pdh*(2 + (-1 + d)*(1 - f_0 - f_1)*mh_2)))) + (-1 + d)*f_0*md_0*(-2*(-1 + pdh)*(2 + (-1 + d)*phh*(1 - f_0 - f_1)*mh_2) + (-1 + d)*f_1*(-((pdh - phh)*mh_1*(2 + (-1 + d)*phh*(1 - f_0 - f_1)*mh_2)) + (-1 + pdh)*md_1*(-2 - (-1 + d)*phh*(1 - f_0 - f_1)*mh_2 + pdh*(2 + (-1 + d)*(1 - f_0 - f_1)*mh_2))))));

gsl_vector_set(f, 0, rhhtplus1);
gsl_vector_set(f, 1, rhdtplus1);
gsl_vector_set(f, 2, rddtplus1);

}


void reproductive_values(void *params, gsl_vector *f)
{
    // 
double d = ((struct rparams *) params)->d;
double md_0 = ((struct rparams *) params)->md_0;
double mh_1 = ((struct rparams *) params)->mh_1;
double md_1 = ((struct rparams *) params)->md_1;
double mh_2 = ((struct rparams *) params)->mh_2;
double pdhinit = ((struct rparams *) params)->pdhinit;
double phhinit = ((struct rparams *) params)->phhinit;
double f_0 = ((struct rparams *) params)->f_0;
double f_1 = ((struct rparams *) params)->f_1;
double pdh = ((struct rparams *) params)->pdh;
double phh = ((struct rparams *) params)->phh;
double rhh = ((struct rparams *) params)->rhh;
double rhd = ((struct rparams *) params)->rhd;
double rdd = ((struct rparams *) params)->rdd;
double vh_1 = ((struct rparams *) params)->vh_1;
double vh_2 = ((struct rparams *) params)->vh_2;
double vd_0 = ((struct rparams *) params)->vd_0;
double vd_1 = ((struct rparams *) params)->vd_1;



    // 
double vh_1_tplus1 = 2/f_1 + (2*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.))/(d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2) + (2*(1 - f_0 - f_1)*(((d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2)*(f_1*(2*d*(-1 + phh)*f_0*md_0 + d*(-1 + phh)*f_1*mh_1) - 2*f_0*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)) - (f_1*(-2*pdh*md_0 + (d*(-2 + 4*pdh - (-2 + 4*pdh)*f_0 - 2*phh*(1 - f_0 - f_1) - pdh*f_1 - phh*f_1)*md_0)/2. + (d*f_1*mh_1)/2. - (d*pdh*f_1*mh_1)/2.) - 2*f_0*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.))*(f_1*(d*(-1 + phh)*f_1*md_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2) - f_1*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)))*(2*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1)*(d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2) - 2*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.)*(-(f_1*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1)) + f_1*((1 + pdh)*md_1 - (-2 + pdh + phh)*mh_1 - 2*d*(-1 + pdh)*f_0*mh_1 + d*(-2 + pdh + phh - 2*(-1 + phh)*(1 - f_0 - f_1) - (-2 + pdh + phh)*f_1)*mh_1 + d*(-1 + pdh)*((-1 + f_1)*md_1 + 2*(1 - f_0 - f_1)*mh_2)))) - (2*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)*(d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2) - 2*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.)*(f_1*(d*(-1 + phh)*f_1*md_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2) - f_1*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)))*((-2*f_0*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1) + f_1*(2*d*(-1 + pdh)*f_0*md_0 + (-3 + 2*pdh + phh)*mh_1 + 2*d*(-1 + pdh)*f_0*mh_1 + d*(3 - 2*pdh - phh + 2*(-1 + phh)*(1 - f_0 - f_1) + (-3 + 2*pdh + phh)*f_1)*mh_1))*(d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2) - (f_1*(-2*pdh*md_0 + (d*(-2 + 4*pdh - (-2 + 4*pdh)*f_0 - 2*phh*(1 - f_0 - f_1) - pdh*f_1 - phh*f_1)*md_0)/2. + (d*f_1*mh_1)/2. - (d*pdh*f_1*mh_1)/2.) - 2*f_0*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.))*(-(f_1*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1)) + f_1*((1 + pdh)*md_1 - (-2 + pdh + phh)*mh_1 - 2*d*(-1 + pdh)*f_0*mh_1 + d*(-2 + pdh + phh - 2*(-1 + phh)*(1 - f_0 - f_1) - (-2 + pdh + phh)*f_1)*mh_1 + d*(-1 + pdh)*((-1 + f_1)*md_1 + 2*(1 - f_0 - f_1)*mh_2))))))/(f_1*(-(((d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2)*(f_1*(-(d*phh*f_1*md_1) - 4*(-1 + phh)*mh_2 - 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-f_0 - f_1)*mh_2 - 4*d*phh*(-f_0 - f_1)*mh_2 - d*(-2 + pdh + phh)*f_1*mh_2) - 2*(1 - f_0 - f_1)*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)) - (-2*(1 - f_0 - f_1)*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.) + f_1*((d*pdh*f_1*md_1)/2. + d*pdh*(1 - f_0 - f_1)*mh_2))*(f_1*(d*(-1 + phh)*f_1*md_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2) - f_1*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)))*((-2*f_0*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1) + f_1*(2*d*(-1 + pdh)*f_0*md_0 + (-3 + 2*pdh + phh)*mh_1 + 2*d*(-1 + pdh)*f_0*mh_1 + d*(3 - 2*pdh - phh + 2*(-1 + phh)*(1 - f_0 - f_1) + (-3 + 2*pdh + phh)*f_1)*mh_1))*(d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2) - (f_1*(-2*pdh*md_0 + (d*(-2 + 4*pdh - (-2 + 4*pdh)*f_0 - 2*phh*(1 - f_0 - f_1) - pdh*f_1 - phh*f_1)*md_0)/2. + (d*f_1*mh_1)/2. - (d*pdh*f_1*mh_1)/2.) - 2*f_0*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.))*(-(f_1*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1)) + f_1*((1 + pdh)*md_1 - (-2 + pdh + phh)*mh_1 - 2*d*(-1 + pdh)*f_0*mh_1 + d*(-2 + pdh + phh - 2*(-1 + phh)*(1 - f_0 - f_1) - (-2 + pdh + phh)*f_1)*mh_1 + d*(-1 + pdh)*((-1 + f_1)*md_1 + 2*(1 - f_0 - f_1)*mh_2))))) + ((d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2)*(f_1*(2*d*(-1 + phh)*f_0*md_0 + d*(-1 + phh)*f_1*mh_1) - 2*f_0*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)) - (f_1*(-2*pdh*md_0 + (d*(-2 + 4*pdh - (-2 + 4*pdh)*f_0 - 2*phh*(1 - f_0 - f_1) - pdh*f_1 - phh*f_1)*md_0)/2. + (d*f_1*mh_1)/2. - (d*pdh*f_1*mh_1)/2.) - 2*f_0*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.))*(f_1*(d*(-1 + phh)*f_1*md_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2) - f_1*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)))*(-((-2*(1 - f_0 - f_1)*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.) + f_1*((d*pdh*f_1*md_1)/2. + d*pdh*(1 - f_0 - f_1)*mh_2))*(-(f_1*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1)) + f_1*((1 + pdh)*md_1 - (-2 + pdh + phh)*mh_1 - 2*d*(-1 + pdh)*f_0*mh_1 + d*(-2 + pdh + phh - 2*(-1 + phh)*(1 - f_0 - f_1) - (-2 + pdh + phh)*f_1)*mh_1 + d*(-1 + pdh)*((-1 + f_1)*md_1 + 2*(1 - f_0 - f_1)*mh_2)))) + (d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2)*(-2*(1 - f_0 - f_1)*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1) + f_1*(-(pdh*md_1) - d*pdh*((-1 + f_1)*md_1 + 2*(1 - f_0 - f_1)*mh_2)))))) - ((-2*(1 - f_0 - f_1)*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.) + f_1*((d*pdh*f_1*md_1)/2. + d*pdh*(1 - f_0 - f_1)*mh_2))*(((d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2)*(f_1*(2*d*(-1 + phh)*f_0*md_0 + d*(-1 + phh)*f_1*mh_1) - 2*f_0*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)) - (f_1*(-2*pdh*md_0 + (d*(-2 + 4*pdh - (-2 + 4*pdh)*f_0 - 2*phh*(1 - f_0 - f_1) - pdh*f_1 - phh*f_1)*md_0)/2. + (d*f_1*mh_1)/2. - (d*pdh*f_1*mh_1)/2.) - 2*f_0*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.))*(f_1*(d*(-1 + phh)*f_1*md_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2) - f_1*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)))*(2*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1)*(d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2) - 2*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.)*(-(f_1*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1)) + f_1*((1 + pdh)*md_1 - (-2 + pdh + phh)*mh_1 - 2*d*(-1 + pdh)*f_0*mh_1 + d*(-2 + pdh + phh - 2*(-1 + phh)*(1 - f_0 - f_1) - (-2 + pdh + phh)*f_1)*mh_1 + d*(-1 + pdh)*((-1 + f_1)*md_1 + 2*(1 - f_0 - f_1)*mh_2)))) - (2*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)*(d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2) - 2*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.)*(f_1*(d*(-1 + phh)*f_1*md_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2) - f_1*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)))*((-2*f_0*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1) + f_1*(2*d*(-1 + pdh)*f_0*md_0 + (-3 + 2*pdh + phh)*mh_1 + 2*d*(-1 + pdh)*f_0*mh_1 + d*(3 - 2*pdh - phh + 2*(-1 + phh)*(1 - f_0 - f_1) + (-3 + 2*pdh + phh)*f_1)*mh_1))*(d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2) - (f_1*(-2*pdh*md_0 + (d*(-2 + 4*pdh - (-2 + 4*pdh)*f_0 - 2*phh*(1 - f_0 - f_1) - pdh*f_1 - phh*f_1)*md_0)/2. + (d*f_1*mh_1)/2. - (d*pdh*f_1*mh_1)/2.) - 2*f_0*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.))*(-(f_1*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1)) + f_1*((1 + pdh)*md_1 - (-2 + pdh + phh)*mh_1 - 2*d*(-1 + pdh)*f_0*mh_1 + d*(-2 + pdh + phh - 2*(-1 + phh)*(1 - f_0 - f_1) - (-2 + pdh + phh)*f_1)*mh_1 + d*(-1 + pdh)*((-1 + f_1)*md_1 + 2*(1 - f_0 - f_1)*mh_2))))))/((d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2)*(-(((d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2)*(f_1*(-(d*phh*f_1*md_1) - 4*(-1 + phh)*mh_2 - 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-f_0 - f_1)*mh_2 - 4*d*phh*(-f_0 - f_1)*mh_2 - d*(-2 + pdh + phh)*f_1*mh_2) - 2*(1 - f_0 - f_1)*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)) - (-2*(1 - f_0 - f_1)*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.) + f_1*((d*pdh*f_1*md_1)/2. + d*pdh*(1 - f_0 - f_1)*mh_2))*(f_1*(d*(-1 + phh)*f_1*md_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2) - f_1*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)))*((-2*f_0*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1) + f_1*(2*d*(-1 + pdh)*f_0*md_0 + (-3 + 2*pdh + phh)*mh_1 + 2*d*(-1 + pdh)*f_0*mh_1 + d*(3 - 2*pdh - phh + 2*(-1 + phh)*(1 - f_0 - f_1) + (-3 + 2*pdh + phh)*f_1)*mh_1))*(d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2) - (f_1*(-2*pdh*md_0 + (d*(-2 + 4*pdh - (-2 + 4*pdh)*f_0 - 2*phh*(1 - f_0 - f_1) - pdh*f_1 - phh*f_1)*md_0)/2. + (d*f_1*mh_1)/2. - (d*pdh*f_1*mh_1)/2.) - 2*f_0*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.))*(-(f_1*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1)) + f_1*((1 + pdh)*md_1 - (-2 + pdh + phh)*mh_1 - 2*d*(-1 + pdh)*f_0*mh_1 + d*(-2 + pdh + phh - 2*(-1 + phh)*(1 - f_0 - f_1) - (-2 + pdh + phh)*f_1)*mh_1 + d*(-1 + pdh)*((-1 + f_1)*md_1 + 2*(1 - f_0 - f_1)*mh_2))))) + ((d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2)*(f_1*(2*d*(-1 + phh)*f_0*md_0 + d*(-1 + phh)*f_1*mh_1) - 2*f_0*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)) - (f_1*(-2*pdh*md_0 + (d*(-2 + 4*pdh - (-2 + 4*pdh)*f_0 - 2*phh*(1 - f_0 - f_1) - pdh*f_1 - phh*f_1)*md_0)/2. + (d*f_1*mh_1)/2. - (d*pdh*f_1*mh_1)/2.) - 2*f_0*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.))*(f_1*(d*(-1 + phh)*f_1*md_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2) - f_1*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)))*(-((-2*(1 - f_0 - f_1)*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.) + f_1*((d*pdh*f_1*md_1)/2. + d*pdh*(1 - f_0 - f_1)*mh_2))*(-(f_1*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1)) + f_1*((1 + pdh)*md_1 - (-2 + pdh + phh)*mh_1 - 2*d*(-1 + pdh)*f_0*mh_1 + d*(-2 + pdh + phh - 2*(-1 + phh)*(1 - f_0 - f_1) - (-2 + pdh + phh)*f_1)*mh_1 + d*(-1 + pdh)*((-1 + f_1)*md_1 + 2*(1 - f_0 - f_1)*mh_2)))) + (d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2)*(-2*(1 - f_0 - f_1)*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1) + f_1*(-(pdh*md_1) - d*pdh*((-1 + f_1)*md_1 + 2*(1 - f_0 - f_1)*mh_2)))))) + (2*f_0*((2*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)*(d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2) - 2*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.)*(f_1*(d*(-1 + phh)*f_1*md_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2) - f_1*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)))/((d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2)*(f_1*(2*d*(-1 + phh)*f_0*md_0 + d*(-1 + phh)*f_1*mh_1) - 2*f_0*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)) - (f_1*(-2*pdh*md_0 + (d*(-2 + 4*pdh - (-2 + 4*pdh)*f_0 - 2*phh*(1 - f_0 - f_1) - pdh*f_1 - phh*f_1)*md_0)/2. + (d*f_1*mh_1)/2. - (d*pdh*f_1*mh_1)/2.) - 2*f_0*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.))*(f_1*(d*(-1 + phh)*f_1*md_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2) - f_1*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2))) - (((d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2)*(f_1*(-(d*phh*f_1*md_1) - 4*(-1 + phh)*mh_2 - 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-f_0 - f_1)*mh_2 - 4*d*phh*(-f_0 - f_1)*mh_2 - d*(-2 + pdh + phh)*f_1*mh_2) - 2*(1 - f_0 - f_1)*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)) - (-2*(1 - f_0 - f_1)*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.) + f_1*((d*pdh*f_1*md_1)/2. + d*pdh*(1 - f_0 - f_1)*mh_2))*(f_1*(d*(-1 + phh)*f_1*md_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2) - f_1*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)))*(((d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2)*(f_1*(2*d*(-1 + phh)*f_0*md_0 + d*(-1 + phh)*f_1*mh_1) - 2*f_0*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)) - (f_1*(-2*pdh*md_0 + (d*(-2 + 4*pdh - (-2 + 4*pdh)*f_0 - 2*phh*(1 - f_0 - f_1) - pdh*f_1 - phh*f_1)*md_0)/2. + (d*f_1*mh_1)/2. - (d*pdh*f_1*mh_1)/2.) - 2*f_0*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.))*(f_1*(d*(-1 + phh)*f_1*md_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2) - f_1*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)))*(2*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1)*(d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2) - 2*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.)*(-(f_1*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1)) + f_1*((1 + pdh)*md_1 - (-2 + pdh + phh)*mh_1 - 2*d*(-1 + pdh)*f_0*mh_1 + d*(-2 + pdh + phh - 2*(-1 + phh)*(1 - f_0 - f_1) - (-2 + pdh + phh)*f_1)*mh_1 + d*(-1 + pdh)*((-1 + f_1)*md_1 + 2*(1 - f_0 - f_1)*mh_2)))) - (2*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)*(d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2) - 2*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.)*(f_1*(d*(-1 + phh)*f_1*md_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2) - f_1*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)))*((-2*f_0*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1) + f_1*(2*d*(-1 + pdh)*f_0*md_0 + (-3 + 2*pdh + phh)*mh_1 + 2*d*(-1 + pdh)*f_0*mh_1 + d*(3 - 2*pdh - phh + 2*(-1 + phh)*(1 - f_0 - f_1) + (-3 + 2*pdh + phh)*f_1)*mh_1))*(d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2) - (f_1*(-2*pdh*md_0 + (d*(-2 + 4*pdh - (-2 + 4*pdh)*f_0 - 2*phh*(1 - f_0 - f_1) - pdh*f_1 - phh*f_1)*md_0)/2. + (d*f_1*mh_1)/2. - (d*pdh*f_1*mh_1)/2.) - 2*f_0*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.))*(-(f_1*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1)) + f_1*((1 + pdh)*md_1 - (-2 + pdh + phh)*mh_1 - 2*d*(-1 + pdh)*f_0*mh_1 + d*(-2 + pdh + phh - 2*(-1 + phh)*(1 - f_0 - f_1) - (-2 + pdh + phh)*f_1)*mh_1 + d*(-1 + pdh)*((-1 + f_1)*md_1 + 2*(1 - f_0 - f_1)*mh_2))))))/(((d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2)*(f_1*(2*d*(-1 + phh)*f_0*md_0 + d*(-1 + phh)*f_1*mh_1) - 2*f_0*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)) - (f_1*(-2*pdh*md_0 + (d*(-2 + 4*pdh - (-2 + 4*pdh)*f_0 - 2*phh*(1 - f_0 - f_1) - pdh*f_1 - phh*f_1)*md_0)/2. + (d*f_1*mh_1)/2. - (d*pdh*f_1*mh_1)/2.) - 2*f_0*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.))*(f_1*(d*(-1 + phh)*f_1*md_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2) - f_1*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)))*(-(((d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2)*(f_1*(-(d*phh*f_1*md_1) - 4*(-1 + phh)*mh_2 - 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-f_0 - f_1)*mh_2 - 4*d*phh*(-f_0 - f_1)*mh_2 - d*(-2 + pdh + phh)*f_1*mh_2) - 2*(1 - f_0 - f_1)*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)) - (-2*(1 - f_0 - f_1)*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.) + f_1*((d*pdh*f_1*md_1)/2. + d*pdh*(1 - f_0 - f_1)*mh_2))*(f_1*(d*(-1 + phh)*f_1*md_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2) - f_1*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)))*((-2*f_0*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1) + f_1*(2*d*(-1 + pdh)*f_0*md_0 + (-3 + 2*pdh + phh)*mh_1 + 2*d*(-1 + pdh)*f_0*mh_1 + d*(3 - 2*pdh - phh + 2*(-1 + phh)*(1 - f_0 - f_1) + (-3 + 2*pdh + phh)*f_1)*mh_1))*(d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2) - (f_1*(-2*pdh*md_0 + (d*(-2 + 4*pdh - (-2 + 4*pdh)*f_0 - 2*phh*(1 - f_0 - f_1) - pdh*f_1 - phh*f_1)*md_0)/2. + (d*f_1*mh_1)/2. - (d*pdh*f_1*mh_1)/2.) - 2*f_0*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.))*(-(f_1*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1)) + f_1*((1 + pdh)*md_1 - (-2 + pdh + phh)*mh_1 - 2*d*(-1 + pdh)*f_0*mh_1 + d*(-2 + pdh + phh - 2*(-1 + phh)*(1 - f_0 - f_1) - (-2 + pdh + phh)*f_1)*mh_1 + d*(-1 + pdh)*((-1 + f_1)*md_1 + 2*(1 - f_0 - f_1)*mh_2))))) + ((d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2)*(f_1*(2*d*(-1 + phh)*f_0*md_0 + d*(-1 + phh)*f_1*mh_1) - 2*f_0*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)) - (f_1*(-2*pdh*md_0 + (d*(-2 + 4*pdh - (-2 + 4*pdh)*f_0 - 2*phh*(1 - f_0 - f_1) - pdh*f_1 - phh*f_1)*md_0)/2. + (d*f_1*mh_1)/2. - (d*pdh*f_1*mh_1)/2.) - 2*f_0*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.))*(f_1*(d*(-1 + phh)*f_1*md_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2) - f_1*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)))*(-((-2*(1 - f_0 - f_1)*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.) + f_1*((d*pdh*f_1*md_1)/2. + d*pdh*(1 - f_0 - f_1)*mh_2))*(-(f_1*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1)) + f_1*((1 + pdh)*md_1 - (-2 + pdh + phh)*mh_1 - 2*d*(-1 + pdh)*f_0*mh_1 + d*(-2 + pdh + phh - 2*(-1 + phh)*(1 - f_0 - f_1) - (-2 + pdh + phh)*f_1)*mh_1 + d*(-1 + pdh)*((-1 + f_1)*md_1 + 2*(1 - f_0 - f_1)*mh_2)))) + (d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2)*(-2*(1 - f_0 - f_1)*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1) + f_1*(-(pdh*md_1) - d*pdh*((-1 + f_1)*md_1 + 2*(1 - f_0 - f_1)*mh_2))))))))/f_1 - ((f_1*(-2*pdh*md_0 + (d*(-2 + 4*pdh - (-2 + 4*pdh)*f_0 - 2*phh*(1 - f_0 - f_1) - pdh*f_1 - phh*f_1)*md_0)/2. + (d*f_1*mh_1)/2. - (d*pdh*f_1*mh_1)/2.) - 2*f_0*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.))*((2*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)*(d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2) - 2*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.)*(f_1*(d*(-1 + phh)*f_1*md_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2) - f_1*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)))/((d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2)*(f_1*(2*d*(-1 + phh)*f_0*md_0 + d*(-1 + phh)*f_1*mh_1) - 2*f_0*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)) - (f_1*(-2*pdh*md_0 + (d*(-2 + 4*pdh - (-2 + 4*pdh)*f_0 - 2*phh*(1 - f_0 - f_1) - pdh*f_1 - phh*f_1)*md_0)/2. + (d*f_1*mh_1)/2. - (d*pdh*f_1*mh_1)/2.) - 2*f_0*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.))*(f_1*(d*(-1 + phh)*f_1*md_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2) - f_1*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2))) - (((d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2)*(f_1*(-(d*phh*f_1*md_1) - 4*(-1 + phh)*mh_2 - 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-f_0 - f_1)*mh_2 - 4*d*phh*(-f_0 - f_1)*mh_2 - d*(-2 + pdh + phh)*f_1*mh_2) - 2*(1 - f_0 - f_1)*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)) - (-2*(1 - f_0 - f_1)*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.) + f_1*((d*pdh*f_1*md_1)/2. + d*pdh*(1 - f_0 - f_1)*mh_2))*(f_1*(d*(-1 + phh)*f_1*md_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2) - f_1*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)))*(((d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2)*(f_1*(2*d*(-1 + phh)*f_0*md_0 + d*(-1 + phh)*f_1*mh_1) - 2*f_0*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)) - (f_1*(-2*pdh*md_0 + (d*(-2 + 4*pdh - (-2 + 4*pdh)*f_0 - 2*phh*(1 - f_0 - f_1) - pdh*f_1 - phh*f_1)*md_0)/2. + (d*f_1*mh_1)/2. - (d*pdh*f_1*mh_1)/2.) - 2*f_0*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.))*(f_1*(d*(-1 + phh)*f_1*md_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2) - f_1*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)))*(2*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1)*(d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2) - 2*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.)*(-(f_1*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1)) + f_1*((1 + pdh)*md_1 - (-2 + pdh + phh)*mh_1 - 2*d*(-1 + pdh)*f_0*mh_1 + d*(-2 + pdh + phh - 2*(-1 + phh)*(1 - f_0 - f_1) - (-2 + pdh + phh)*f_1)*mh_1 + d*(-1 + pdh)*((-1 + f_1)*md_1 + 2*(1 - f_0 - f_1)*mh_2)))) - (2*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)*(d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2) - 2*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.)*(f_1*(d*(-1 + phh)*f_1*md_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2) - f_1*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)))*((-2*f_0*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1) + f_1*(2*d*(-1 + pdh)*f_0*md_0 + (-3 + 2*pdh + phh)*mh_1 + 2*d*(-1 + pdh)*f_0*mh_1 + d*(3 - 2*pdh - phh + 2*(-1 + phh)*(1 - f_0 - f_1) + (-3 + 2*pdh + phh)*f_1)*mh_1))*(d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2) - (f_1*(-2*pdh*md_0 + (d*(-2 + 4*pdh - (-2 + 4*pdh)*f_0 - 2*phh*(1 - f_0 - f_1) - pdh*f_1 - phh*f_1)*md_0)/2. + (d*f_1*mh_1)/2. - (d*pdh*f_1*mh_1)/2.) - 2*f_0*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.))*(-(f_1*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1)) + f_1*((1 + pdh)*md_1 - (-2 + pdh + phh)*mh_1 - 2*d*(-1 + pdh)*f_0*mh_1 + d*(-2 + pdh + phh - 2*(-1 + phh)*(1 - f_0 - f_1) - (-2 + pdh + phh)*f_1)*mh_1 + d*(-1 + pdh)*((-1 + f_1)*md_1 + 2*(1 - f_0 - f_1)*mh_2))))))/(((d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2)*(f_1*(2*d*(-1 + phh)*f_0*md_0 + d*(-1 + phh)*f_1*mh_1) - 2*f_0*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)) - (f_1*(-2*pdh*md_0 + (d*(-2 + 4*pdh - (-2 + 4*pdh)*f_0 - 2*phh*(1 - f_0 - f_1) - pdh*f_1 - phh*f_1)*md_0)/2. + (d*f_1*mh_1)/2. - (d*pdh*f_1*mh_1)/2.) - 2*f_0*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.))*(f_1*(d*(-1 + phh)*f_1*md_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2) - f_1*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)))*(-(((d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2)*(f_1*(-(d*phh*f_1*md_1) - 4*(-1 + phh)*mh_2 - 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-f_0 - f_1)*mh_2 - 4*d*phh*(-f_0 - f_1)*mh_2 - d*(-2 + pdh + phh)*f_1*mh_2) - 2*(1 - f_0 - f_1)*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)) - (-2*(1 - f_0 - f_1)*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.) + f_1*((d*pdh*f_1*md_1)/2. + d*pdh*(1 - f_0 - f_1)*mh_2))*(f_1*(d*(-1 + phh)*f_1*md_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2) - f_1*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)))*((-2*f_0*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1) + f_1*(2*d*(-1 + pdh)*f_0*md_0 + (-3 + 2*pdh + phh)*mh_1 + 2*d*(-1 + pdh)*f_0*mh_1 + d*(3 - 2*pdh - phh + 2*(-1 + phh)*(1 - f_0 - f_1) + (-3 + 2*pdh + phh)*f_1)*mh_1))*(d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2) - (f_1*(-2*pdh*md_0 + (d*(-2 + 4*pdh - (-2 + 4*pdh)*f_0 - 2*phh*(1 - f_0 - f_1) - pdh*f_1 - phh*f_1)*md_0)/2. + (d*f_1*mh_1)/2. - (d*pdh*f_1*mh_1)/2.) - 2*f_0*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.))*(-(f_1*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1)) + f_1*((1 + pdh)*md_1 - (-2 + pdh + phh)*mh_1 - 2*d*(-1 + pdh)*f_0*mh_1 + d*(-2 + pdh + phh - 2*(-1 + phh)*(1 - f_0 - f_1) - (-2 + pdh + phh)*f_1)*mh_1 + d*(-1 + pdh)*((-1 + f_1)*md_1 + 2*(1 - f_0 - f_1)*mh_2))))) + ((d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2)*(f_1*(2*d*(-1 + phh)*f_0*md_0 + d*(-1 + phh)*f_1*mh_1) - 2*f_0*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)) - (f_1*(-2*pdh*md_0 + (d*(-2 + 4*pdh - (-2 + 4*pdh)*f_0 - 2*phh*(1 - f_0 - f_1) - pdh*f_1 - phh*f_1)*md_0)/2. + (d*f_1*mh_1)/2. - (d*pdh*f_1*mh_1)/2.) - 2*f_0*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.))*(f_1*(d*(-1 + phh)*f_1*md_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2) - f_1*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)))*(-((-2*(1 - f_0 - f_1)*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.) + f_1*((d*pdh*f_1*md_1)/2. + d*pdh*(1 - f_0 - f_1)*mh_2))*(-(f_1*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1)) + f_1*((1 + pdh)*md_1 - (-2 + pdh + phh)*mh_1 - 2*d*(-1 + pdh)*f_0*mh_1 + d*(-2 + pdh + phh - 2*(-1 + phh)*(1 - f_0 - f_1) - (-2 + pdh + phh)*f_1)*mh_1 + d*(-1 + pdh)*((-1 + f_1)*md_1 + 2*(1 - f_0 - f_1)*mh_2)))) + (d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2)*(-2*(1 - f_0 - f_1)*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1) + f_1*(-(pdh*md_1) - d*pdh*((-1 + f_1)*md_1 + 2*(1 - f_0 - f_1)*mh_2))))))))/(d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2);

double vh_2_tplus1 = -((((d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2)*(f_1*(2*d*(-1 + phh)*f_0*md_0 + d*(-1 + phh)*f_1*mh_1) - 2*f_0*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)) - (f_1*(-2*pdh*md_0 + (d*(-2 + 4*pdh - (-2 + 4*pdh)*f_0 - 2*phh*(1 - f_0 - f_1) - pdh*f_1 - phh*f_1)*md_0)/2. + (d*f_1*mh_1)/2. - (d*pdh*f_1*mh_1)/2.) - 2*f_0*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.))*(f_1*(d*(-1 + phh)*f_1*md_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2) - f_1*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)))*(2*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1)*(d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2) - 2*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.)*(-(f_1*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1)) + f_1*((1 + pdh)*md_1 - (-2 + pdh + phh)*mh_1 - 2*d*(-1 + pdh)*f_0*mh_1 + d*(-2 + pdh + phh - 2*(-1 + phh)*(1 - f_0 - f_1) - (-2 + pdh + phh)*f_1)*mh_1 + d*(-1 + pdh)*((-1 + f_1)*md_1 + 2*(1 - f_0 - f_1)*mh_2)))) - (2*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)*(d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2) - 2*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.)*(f_1*(d*(-1 + phh)*f_1*md_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2) - f_1*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)))*((-2*f_0*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1) + f_1*(2*d*(-1 + pdh)*f_0*md_0 + (-3 + 2*pdh + phh)*mh_1 + 2*d*(-1 + pdh)*f_0*mh_1 + d*(3 - 2*pdh - phh + 2*(-1 + phh)*(1 - f_0 - f_1) + (-3 + 2*pdh + phh)*f_1)*mh_1))*(d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2) - (f_1*(-2*pdh*md_0 + (d*(-2 + 4*pdh - (-2 + 4*pdh)*f_0 - 2*phh*(1 - f_0 - f_1) - pdh*f_1 - phh*f_1)*md_0)/2. + (d*f_1*mh_1)/2. - (d*pdh*f_1*mh_1)/2.) - 2*f_0*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.))*(-(f_1*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1)) + f_1*((1 + pdh)*md_1 - (-2 + pdh + phh)*mh_1 - 2*d*(-1 + pdh)*f_0*mh_1 + d*(-2 + pdh + phh - 2*(-1 + phh)*(1 - f_0 - f_1) - (-2 + pdh + phh)*f_1)*mh_1 + d*(-1 + pdh)*((-1 + f_1)*md_1 + 2*(1 - f_0 - f_1)*mh_2)))))/(-(((d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2)*(f_1*(-(d*phh*f_1*md_1) - 4*(-1 + phh)*mh_2 - 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-f_0 - f_1)*mh_2 - 4*d*phh*(-f_0 - f_1)*mh_2 - d*(-2 + pdh + phh)*f_1*mh_2) - 2*(1 - f_0 - f_1)*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)) - (-2*(1 - f_0 - f_1)*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.) + f_1*((d*pdh*f_1*md_1)/2. + d*pdh*(1 - f_0 - f_1)*mh_2))*(f_1*(d*(-1 + phh)*f_1*md_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2) - f_1*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)))*((-2*f_0*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1) + f_1*(2*d*(-1 + pdh)*f_0*md_0 + (-3 + 2*pdh + phh)*mh_1 + 2*d*(-1 + pdh)*f_0*mh_1 + d*(3 - 2*pdh - phh + 2*(-1 + phh)*(1 - f_0 - f_1) + (-3 + 2*pdh + phh)*f_1)*mh_1))*(d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2) - (f_1*(-2*pdh*md_0 + (d*(-2 + 4*pdh - (-2 + 4*pdh)*f_0 - 2*phh*(1 - f_0 - f_1) - pdh*f_1 - phh*f_1)*md_0)/2. + (d*f_1*mh_1)/2. - (d*pdh*f_1*mh_1)/2.) - 2*f_0*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.))*(-(f_1*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1)) + f_1*((1 + pdh)*md_1 - (-2 + pdh + phh)*mh_1 - 2*d*(-1 + pdh)*f_0*mh_1 + d*(-2 + pdh + phh - 2*(-1 + phh)*(1 - f_0 - f_1) - (-2 + pdh + phh)*f_1)*mh_1 + d*(-1 + pdh)*((-1 + f_1)*md_1 + 2*(1 - f_0 - f_1)*mh_2))))) + ((d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2)*(f_1*(2*d*(-1 + phh)*f_0*md_0 + d*(-1 + phh)*f_1*mh_1) - 2*f_0*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)) - (f_1*(-2*pdh*md_0 + (d*(-2 + 4*pdh - (-2 + 4*pdh)*f_0 - 2*phh*(1 - f_0 - f_1) - pdh*f_1 - phh*f_1)*md_0)/2. + (d*f_1*mh_1)/2. - (d*pdh*f_1*mh_1)/2.) - 2*f_0*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.))*(f_1*(d*(-1 + phh)*f_1*md_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2) - f_1*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)))*(-((-2*(1 - f_0 - f_1)*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.) + f_1*((d*pdh*f_1*md_1)/2. + d*pdh*(1 - f_0 - f_1)*mh_2))*(-(f_1*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1)) + f_1*((1 + pdh)*md_1 - (-2 + pdh + phh)*mh_1 - 2*d*(-1 + pdh)*f_0*mh_1 + d*(-2 + pdh + phh - 2*(-1 + phh)*(1 - f_0 - f_1) - (-2 + pdh + phh)*f_1)*mh_1 + d*(-1 + pdh)*((-1 + f_1)*md_1 + 2*(1 - f_0 - f_1)*mh_2)))) + (d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2)*(-2*(1 - f_0 - f_1)*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1) + f_1*(-(pdh*md_1) - d*pdh*((-1 + f_1)*md_1 + 2*(1 - f_0 - f_1)*mh_2))))));

double vd_0_tplus1 = -((2*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)*(d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2) - 2*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.)*(f_1*(d*(-1 + phh)*f_1*md_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2) - f_1*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)))/((d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2)*(f_1*(2*d*(-1 + phh)*f_0*md_0 + d*(-1 + phh)*f_1*mh_1) - 2*f_0*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)) - (f_1*(-2*pdh*md_0 + (d*(-2 + 4*pdh - (-2 + 4*pdh)*f_0 - 2*phh*(1 - f_0 - f_1) - pdh*f_1 - phh*f_1)*md_0)/2. + (d*f_1*mh_1)/2. - (d*pdh*f_1*mh_1)/2.) - 2*f_0*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.))*(f_1*(d*(-1 + phh)*f_1*md_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2) - f_1*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)))) + (((d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2)*(f_1*(-(d*phh*f_1*md_1) - 4*(-1 + phh)*mh_2 - 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-f_0 - f_1)*mh_2 - 4*d*phh*(-f_0 - f_1)*mh_2 - d*(-2 + pdh + phh)*f_1*mh_2) - 2*(1 - f_0 - f_1)*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)) - (-2*(1 - f_0 - f_1)*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.) + f_1*((d*pdh*f_1*md_1)/2. + d*pdh*(1 - f_0 - f_1)*mh_2))*(f_1*(d*(-1 + phh)*f_1*md_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2) - f_1*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)))*(((d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2)*(f_1*(2*d*(-1 + phh)*f_0*md_0 + d*(-1 + phh)*f_1*mh_1) - 2*f_0*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)) - (f_1*(-2*pdh*md_0 + (d*(-2 + 4*pdh - (-2 + 4*pdh)*f_0 - 2*phh*(1 - f_0 - f_1) - pdh*f_1 - phh*f_1)*md_0)/2. + (d*f_1*mh_1)/2. - (d*pdh*f_1*mh_1)/2.) - 2*f_0*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.))*(f_1*(d*(-1 + phh)*f_1*md_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2) - f_1*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)))*(2*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1)*(d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2) - 2*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.)*(-(f_1*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1)) + f_1*((1 + pdh)*md_1 - (-2 + pdh + phh)*mh_1 - 2*d*(-1 + pdh)*f_0*mh_1 + d*(-2 + pdh + phh - 2*(-1 + phh)*(1 - f_0 - f_1) - (-2 + pdh + phh)*f_1)*mh_1 + d*(-1 + pdh)*((-1 + f_1)*md_1 + 2*(1 - f_0 - f_1)*mh_2)))) - (2*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)*(d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2) - 2*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.)*(f_1*(d*(-1 + phh)*f_1*md_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2) - f_1*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)))*((-2*f_0*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1) + f_1*(2*d*(-1 + pdh)*f_0*md_0 + (-3 + 2*pdh + phh)*mh_1 + 2*d*(-1 + pdh)*f_0*mh_1 + d*(3 - 2*pdh - phh + 2*(-1 + phh)*(1 - f_0 - f_1) + (-3 + 2*pdh + phh)*f_1)*mh_1))*(d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2) - (f_1*(-2*pdh*md_0 + (d*(-2 + 4*pdh - (-2 + 4*pdh)*f_0 - 2*phh*(1 - f_0 - f_1) - pdh*f_1 - phh*f_1)*md_0)/2. + (d*f_1*mh_1)/2. - (d*pdh*f_1*mh_1)/2.) - 2*f_0*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.))*(-(f_1*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1)) + f_1*((1 + pdh)*md_1 - (-2 + pdh + phh)*mh_1 - 2*d*(-1 + pdh)*f_0*mh_1 + d*(-2 + pdh + phh - 2*(-1 + phh)*(1 - f_0 - f_1) - (-2 + pdh + phh)*f_1)*mh_1 + d*(-1 + pdh)*((-1 + f_1)*md_1 + 2*(1 - f_0 - f_1)*mh_2))))))/(((d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2)*(f_1*(2*d*(-1 + phh)*f_0*md_0 + d*(-1 + phh)*f_1*mh_1) - 2*f_0*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)) - (f_1*(-2*pdh*md_0 + (d*(-2 + 4*pdh - (-2 + 4*pdh)*f_0 - 2*phh*(1 - f_0 - f_1) - pdh*f_1 - phh*f_1)*md_0)/2. + (d*f_1*mh_1)/2. - (d*pdh*f_1*mh_1)/2.) - 2*f_0*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.))*(f_1*(d*(-1 + phh)*f_1*md_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2) - f_1*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)))*(-(((d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2)*(f_1*(-(d*phh*f_1*md_1) - 4*(-1 + phh)*mh_2 - 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-f_0 - f_1)*mh_2 - 4*d*phh*(-f_0 - f_1)*mh_2 - d*(-2 + pdh + phh)*f_1*mh_2) - 2*(1 - f_0 - f_1)*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)) - (-2*(1 - f_0 - f_1)*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.) + f_1*((d*pdh*f_1*md_1)/2. + d*pdh*(1 - f_0 - f_1)*mh_2))*(f_1*(d*(-1 + phh)*f_1*md_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2) - f_1*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)))*((-2*f_0*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1) + f_1*(2*d*(-1 + pdh)*f_0*md_0 + (-3 + 2*pdh + phh)*mh_1 + 2*d*(-1 + pdh)*f_0*mh_1 + d*(3 - 2*pdh - phh + 2*(-1 + phh)*(1 - f_0 - f_1) + (-3 + 2*pdh + phh)*f_1)*mh_1))*(d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2) - (f_1*(-2*pdh*md_0 + (d*(-2 + 4*pdh - (-2 + 4*pdh)*f_0 - 2*phh*(1 - f_0 - f_1) - pdh*f_1 - phh*f_1)*md_0)/2. + (d*f_1*mh_1)/2. - (d*pdh*f_1*mh_1)/2.) - 2*f_0*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.))*(-(f_1*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1)) + f_1*((1 + pdh)*md_1 - (-2 + pdh + phh)*mh_1 - 2*d*(-1 + pdh)*f_0*mh_1 + d*(-2 + pdh + phh - 2*(-1 + phh)*(1 - f_0 - f_1) - (-2 + pdh + phh)*f_1)*mh_1 + d*(-1 + pdh)*((-1 + f_1)*md_1 + 2*(1 - f_0 - f_1)*mh_2))))) + ((d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2)*(f_1*(2*d*(-1 + phh)*f_0*md_0 + d*(-1 + phh)*f_1*mh_1) - 2*f_0*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)) - (f_1*(-2*pdh*md_0 + (d*(-2 + 4*pdh - (-2 + 4*pdh)*f_0 - 2*phh*(1 - f_0 - f_1) - pdh*f_1 - phh*f_1)*md_0)/2. + (d*f_1*mh_1)/2. - (d*pdh*f_1*mh_1)/2.) - 2*f_0*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.))*(f_1*(d*(-1 + phh)*f_1*md_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2) - f_1*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)))*(-((-2*(1 - f_0 - f_1)*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.) + f_1*((d*pdh*f_1*md_1)/2. + d*pdh*(1 - f_0 - f_1)*mh_2))*(-(f_1*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1)) + f_1*((1 + pdh)*md_1 - (-2 + pdh + phh)*mh_1 - 2*d*(-1 + pdh)*f_0*mh_1 + d*(-2 + pdh + phh - 2*(-1 + phh)*(1 - f_0 - f_1) - (-2 + pdh + phh)*f_1)*mh_1 + d*(-1 + pdh)*((-1 + f_1)*md_1 + 2*(1 - f_0 - f_1)*mh_2)))) + (d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2)*(-2*(1 - f_0 - f_1)*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1) + f_1*(-(pdh*md_1) - d*pdh*((-1 + f_1)*md_1 + 2*(1 - f_0 - f_1)*mh_2))))));

double vd_1_tplus1 = (-2*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.))/(d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2) + ((-2*(1 - f_0 - f_1)*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.) + f_1*((d*pdh*f_1*md_1)/2. + d*pdh*(1 - f_0 - f_1)*mh_2))*(((d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2)*(f_1*(2*d*(-1 + phh)*f_0*md_0 + d*(-1 + phh)*f_1*mh_1) - 2*f_0*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)) - (f_1*(-2*pdh*md_0 + (d*(-2 + 4*pdh - (-2 + 4*pdh)*f_0 - 2*phh*(1 - f_0 - f_1) - pdh*f_1 - phh*f_1)*md_0)/2. + (d*f_1*mh_1)/2. - (d*pdh*f_1*mh_1)/2.) - 2*f_0*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.))*(f_1*(d*(-1 + phh)*f_1*md_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2) - f_1*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)))*(2*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1)*(d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2) - 2*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.)*(-(f_1*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1)) + f_1*((1 + pdh)*md_1 - (-2 + pdh + phh)*mh_1 - 2*d*(-1 + pdh)*f_0*mh_1 + d*(-2 + pdh + phh - 2*(-1 + phh)*(1 - f_0 - f_1) - (-2 + pdh + phh)*f_1)*mh_1 + d*(-1 + pdh)*((-1 + f_1)*md_1 + 2*(1 - f_0 - f_1)*mh_2)))) - (2*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)*(d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2) - 2*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.)*(f_1*(d*(-1 + phh)*f_1*md_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2) - f_1*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)))*((-2*f_0*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1) + f_1*(2*d*(-1 + pdh)*f_0*md_0 + (-3 + 2*pdh + phh)*mh_1 + 2*d*(-1 + pdh)*f_0*mh_1 + d*(3 - 2*pdh - phh + 2*(-1 + phh)*(1 - f_0 - f_1) + (-3 + 2*pdh + phh)*f_1)*mh_1))*(d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2) - (f_1*(-2*pdh*md_0 + (d*(-2 + 4*pdh - (-2 + 4*pdh)*f_0 - 2*phh*(1 - f_0 - f_1) - pdh*f_1 - phh*f_1)*md_0)/2. + (d*f_1*mh_1)/2. - (d*pdh*f_1*mh_1)/2.) - 2*f_0*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.))*(-(f_1*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1)) + f_1*((1 + pdh)*md_1 - (-2 + pdh + phh)*mh_1 - 2*d*(-1 + pdh)*f_0*mh_1 + d*(-2 + pdh + phh - 2*(-1 + phh)*(1 - f_0 - f_1) - (-2 + pdh + phh)*f_1)*mh_1 + d*(-1 + pdh)*((-1 + f_1)*md_1 + 2*(1 - f_0 - f_1)*mh_2))))))/((d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2)*(-(((d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2)*(f_1*(-(d*phh*f_1*md_1) - 4*(-1 + phh)*mh_2 - 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-f_0 - f_1)*mh_2 - 4*d*phh*(-f_0 - f_1)*mh_2 - d*(-2 + pdh + phh)*f_1*mh_2) - 2*(1 - f_0 - f_1)*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)) - (-2*(1 - f_0 - f_1)*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.) + f_1*((d*pdh*f_1*md_1)/2. + d*pdh*(1 - f_0 - f_1)*mh_2))*(f_1*(d*(-1 + phh)*f_1*md_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2) - f_1*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)))*((-2*f_0*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1) + f_1*(2*d*(-1 + pdh)*f_0*md_0 + (-3 + 2*pdh + phh)*mh_1 + 2*d*(-1 + pdh)*f_0*mh_1 + d*(3 - 2*pdh - phh + 2*(-1 + phh)*(1 - f_0 - f_1) + (-3 + 2*pdh + phh)*f_1)*mh_1))*(d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2) - (f_1*(-2*pdh*md_0 + (d*(-2 + 4*pdh - (-2 + 4*pdh)*f_0 - 2*phh*(1 - f_0 - f_1) - pdh*f_1 - phh*f_1)*md_0)/2. + (d*f_1*mh_1)/2. - (d*pdh*f_1*mh_1)/2.) - 2*f_0*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.))*(-(f_1*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1)) + f_1*((1 + pdh)*md_1 - (-2 + pdh + phh)*mh_1 - 2*d*(-1 + pdh)*f_0*mh_1 + d*(-2 + pdh + phh - 2*(-1 + phh)*(1 - f_0 - f_1) - (-2 + pdh + phh)*f_1)*mh_1 + d*(-1 + pdh)*((-1 + f_1)*md_1 + 2*(1 - f_0 - f_1)*mh_2))))) + ((d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2)*(f_1*(2*d*(-1 + phh)*f_0*md_0 + d*(-1 + phh)*f_1*mh_1) - 2*f_0*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)) - (f_1*(-2*pdh*md_0 + (d*(-2 + 4*pdh - (-2 + 4*pdh)*f_0 - 2*phh*(1 - f_0 - f_1) - pdh*f_1 - phh*f_1)*md_0)/2. + (d*f_1*mh_1)/2. - (d*pdh*f_1*mh_1)/2.) - 2*f_0*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.))*(f_1*(d*(-1 + phh)*f_1*md_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2) - f_1*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)))*(-((-2*(1 - f_0 - f_1)*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.) + f_1*((d*pdh*f_1*md_1)/2. + d*pdh*(1 - f_0 - f_1)*mh_2))*(-(f_1*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1)) + f_1*((1 + pdh)*md_1 - (-2 + pdh + phh)*mh_1 - 2*d*(-1 + pdh)*f_0*mh_1 + d*(-2 + pdh + phh - 2*(-1 + phh)*(1 - f_0 - f_1) - (-2 + pdh + phh)*f_1)*mh_1 + d*(-1 + pdh)*((-1 + f_1)*md_1 + 2*(1 - f_0 - f_1)*mh_2)))) + (d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2)*(-2*(1 - f_0 - f_1)*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1) + f_1*(-(pdh*md_1) - d*pdh*((-1 + f_1)*md_1 + 2*(1 - f_0 - f_1)*mh_2)))))) + ((f_1*(-2*pdh*md_0 + (d*(-2 + 4*pdh - (-2 + 4*pdh)*f_0 - 2*phh*(1 - f_0 - f_1) - pdh*f_1 - phh*f_1)*md_0)/2. + (d*f_1*mh_1)/2. - (d*pdh*f_1*mh_1)/2.) - 2*f_0*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.))*((2*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)*(d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2) - 2*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.)*(f_1*(d*(-1 + phh)*f_1*md_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2) - f_1*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)))/((d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2)*(f_1*(2*d*(-1 + phh)*f_0*md_0 + d*(-1 + phh)*f_1*mh_1) - 2*f_0*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)) - (f_1*(-2*pdh*md_0 + (d*(-2 + 4*pdh - (-2 + 4*pdh)*f_0 - 2*phh*(1 - f_0 - f_1) - pdh*f_1 - phh*f_1)*md_0)/2. + (d*f_1*mh_1)/2. - (d*pdh*f_1*mh_1)/2.) - 2*f_0*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.))*(f_1*(d*(-1 + phh)*f_1*md_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2) - f_1*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2))) - (((d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2)*(f_1*(-(d*phh*f_1*md_1) - 4*(-1 + phh)*mh_2 - 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-f_0 - f_1)*mh_2 - 4*d*phh*(-f_0 - f_1)*mh_2 - d*(-2 + pdh + phh)*f_1*mh_2) - 2*(1 - f_0 - f_1)*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)) - (-2*(1 - f_0 - f_1)*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.) + f_1*((d*pdh*f_1*md_1)/2. + d*pdh*(1 - f_0 - f_1)*mh_2))*(f_1*(d*(-1 + phh)*f_1*md_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2) - f_1*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)))*(((d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2)*(f_1*(2*d*(-1 + phh)*f_0*md_0 + d*(-1 + phh)*f_1*mh_1) - 2*f_0*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)) - (f_1*(-2*pdh*md_0 + (d*(-2 + 4*pdh - (-2 + 4*pdh)*f_0 - 2*phh*(1 - f_0 - f_1) - pdh*f_1 - phh*f_1)*md_0)/2. + (d*f_1*mh_1)/2. - (d*pdh*f_1*mh_1)/2.) - 2*f_0*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.))*(f_1*(d*(-1 + phh)*f_1*md_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2) - f_1*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)))*(2*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1)*(d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2) - 2*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.)*(-(f_1*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1)) + f_1*((1 + pdh)*md_1 - (-2 + pdh + phh)*mh_1 - 2*d*(-1 + pdh)*f_0*mh_1 + d*(-2 + pdh + phh - 2*(-1 + phh)*(1 - f_0 - f_1) - (-2 + pdh + phh)*f_1)*mh_1 + d*(-1 + pdh)*((-1 + f_1)*md_1 + 2*(1 - f_0 - f_1)*mh_2)))) - (2*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)*(d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2) - 2*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.)*(f_1*(d*(-1 + phh)*f_1*md_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2) - f_1*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)))*((-2*f_0*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1) + f_1*(2*d*(-1 + pdh)*f_0*md_0 + (-3 + 2*pdh + phh)*mh_1 + 2*d*(-1 + pdh)*f_0*mh_1 + d*(3 - 2*pdh - phh + 2*(-1 + phh)*(1 - f_0 - f_1) + (-3 + 2*pdh + phh)*f_1)*mh_1))*(d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2) - (f_1*(-2*pdh*md_0 + (d*(-2 + 4*pdh - (-2 + 4*pdh)*f_0 - 2*phh*(1 - f_0 - f_1) - pdh*f_1 - phh*f_1)*md_0)/2. + (d*f_1*mh_1)/2. - (d*pdh*f_1*mh_1)/2.) - 2*f_0*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.))*(-(f_1*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1)) + f_1*((1 + pdh)*md_1 - (-2 + pdh + phh)*mh_1 - 2*d*(-1 + pdh)*f_0*mh_1 + d*(-2 + pdh + phh - 2*(-1 + phh)*(1 - f_0 - f_1) - (-2 + pdh + phh)*f_1)*mh_1 + d*(-1 + pdh)*((-1 + f_1)*md_1 + 2*(1 - f_0 - f_1)*mh_2))))))/(((d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2)*(f_1*(2*d*(-1 + phh)*f_0*md_0 + d*(-1 + phh)*f_1*mh_1) - 2*f_0*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)) - (f_1*(-2*pdh*md_0 + (d*(-2 + 4*pdh - (-2 + 4*pdh)*f_0 - 2*phh*(1 - f_0 - f_1) - pdh*f_1 - phh*f_1)*md_0)/2. + (d*f_1*mh_1)/2. - (d*pdh*f_1*mh_1)/2.) - 2*f_0*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.))*(f_1*(d*(-1 + phh)*f_1*md_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2) - f_1*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)))*(-(((d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2)*(f_1*(-(d*phh*f_1*md_1) - 4*(-1 + phh)*mh_2 - 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-f_0 - f_1)*mh_2 - 4*d*phh*(-f_0 - f_1)*mh_2 - d*(-2 + pdh + phh)*f_1*mh_2) - 2*(1 - f_0 - f_1)*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)) - (-2*(1 - f_0 - f_1)*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.) + f_1*((d*pdh*f_1*md_1)/2. + d*pdh*(1 - f_0 - f_1)*mh_2))*(f_1*(d*(-1 + phh)*f_1*md_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2) - f_1*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)))*((-2*f_0*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1) + f_1*(2*d*(-1 + pdh)*f_0*md_0 + (-3 + 2*pdh + phh)*mh_1 + 2*d*(-1 + pdh)*f_0*mh_1 + d*(3 - 2*pdh - phh + 2*(-1 + phh)*(1 - f_0 - f_1) + (-3 + 2*pdh + phh)*f_1)*mh_1))*(d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2) - (f_1*(-2*pdh*md_0 + (d*(-2 + 4*pdh - (-2 + 4*pdh)*f_0 - 2*phh*(1 - f_0 - f_1) - pdh*f_1 - phh*f_1)*md_0)/2. + (d*f_1*mh_1)/2. - (d*pdh*f_1*mh_1)/2.) - 2*f_0*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.))*(-(f_1*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1)) + f_1*((1 + pdh)*md_1 - (-2 + pdh + phh)*mh_1 - 2*d*(-1 + pdh)*f_0*mh_1 + d*(-2 + pdh + phh - 2*(-1 + phh)*(1 - f_0 - f_1) - (-2 + pdh + phh)*f_1)*mh_1 + d*(-1 + pdh)*((-1 + f_1)*md_1 + 2*(1 - f_0 - f_1)*mh_2))))) + ((d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2)*(f_1*(2*d*(-1 + phh)*f_0*md_0 + d*(-1 + phh)*f_1*mh_1) - 2*f_0*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)) - (f_1*(-2*pdh*md_0 + (d*(-2 + 4*pdh - (-2 + 4*pdh)*f_0 - 2*phh*(1 - f_0 - f_1) - pdh*f_1 - phh*f_1)*md_0)/2. + (d*f_1*mh_1)/2. - (d*pdh*f_1*mh_1)/2.) - 2*f_0*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.))*(f_1*(d*(-1 + phh)*f_1*md_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2) - f_1*(-2*d*phh*f_0*md_0 - d*phh*f_1*mh_1 + 2*(-1 + phh)*mh_2 + 2*d*(-1 + pdh)*f_0*mh_2 + 2*d*(-1 + phh)*(-f_0 - f_1)*mh_2 + d*(-2 + pdh + phh)*f_1*mh_2)))*(-((-2*(1 - f_0 - f_1)*(pdh*md_0 + d*pdh*(-1 + f_0)*md_0 + (d*pdh*f_1*mh_1)/2.) + f_1*((d*pdh*f_1*md_1)/2. + d*pdh*(1 - f_0 - f_1)*mh_2))*(-(f_1*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1)) + f_1*((1 + pdh)*md_1 - (-2 + pdh + phh)*mh_1 - 2*d*(-1 + pdh)*f_0*mh_1 + d*(-2 + pdh + phh - 2*(-1 + phh)*(1 - f_0 - f_1) - (-2 + pdh + phh)*f_1)*mh_1 + d*(-1 + pdh)*((-1 + f_1)*md_1 + 2*(1 - f_0 - f_1)*mh_2)))) + (d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2)*(-2*(1 - f_0 - f_1)*(-2*d*pdh*f_0*md_0 - pdh*mh_1 - d*pdh*(-1 + f_1)*mh_1) + f_1*(-(pdh*md_1) - d*pdh*((-1 + f_1)*md_1 + 2*(1 - f_0 - f_1)*mh_2))))))))/(d*phh*(1 - f_0 - f_1)*f_1*md_0 + (d*pdh*pow(f_1,2)*md_0)/2. + (d*phh*pow(f_1,2)*md_0)/2. + (d*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*md_1)/2. - (d*pdh*pow(f_1,2)*mh_1)/2. + d*(1 - f_0 - f_1)*f_1*mh_2 - d*pdh*(1 - f_0 - f_1)*f_1*mh_2);

gsl_vector_set(f, 0, vh_1_tplus1);
gsl_vector_set(f, 1, vh_2_tplus1);
gsl_vector_set(f, 2, vd_0_tplus1);
gsl_vector_set(f, 3, vd_1_tplus1);
}

void selgrads(void *params, gsl_vector *f)
{
    // 
double d = ((struct rparams *) params)->d;
double md_0 = ((struct rparams *) params)->md_0;
double mh_1 = ((struct rparams *) params)->mh_1;
double md_1 = ((struct rparams *) params)->md_1;
double mh_2 = ((struct rparams *) params)->mh_2;
double pdhinit = ((struct rparams *) params)->pdhinit;
double phhinit = ((struct rparams *) params)->phhinit;
double f_0 = ((struct rparams *) params)->f_0;
double f_1 = ((struct rparams *) params)->f_1;
double pdh = ((struct rparams *) params)->pdh;
double phh = ((struct rparams *) params)->phh;
double rhh = ((struct rparams *) params)->rhh;
double rhd = ((struct rparams *) params)->rhd;
double rdd = ((struct rparams *) params)->rdd;
double vh_1 = ((struct rparams *) params)->vh_1;
double vh_2 = ((struct rparams *) params)->vh_2;
double vd_0 = ((struct rparams *) params)->vd_0;
double vd_1 = ((struct rparams *) params)->vd_1;


    // 
double phhtplus1=bound(phh + 0.01*((1 - f_0 - f_1)*(d*f_0*md_0*(-vd_0 + vh_1) + (d*f_1*mh_1*(-vd_0 + vh_1))/2. - ((1 - d)*mh_2*(vd_1 - vh_2))/2. - ((1 - d)*rhh*mh_2*(vh_1 - vh_2))/2. - ((1 - d)*mh_2*(vd_1 + vh_1 - vh_2))/2. + ((1 - d)*mh_2*vh_2)/2. + (d*f_1*md_1*(-vd_1 + vh_2))/2. + d*(1 - f_0 - f_1)*mh_2*(-vd_1 + vh_2)) + (f_1*(-((1 - d)*rhd*mh_1*(vd_0 - vd_1))/2. - ((1 - d)*md_1*vd_1)/2. - ((1 - d)*mh_1*(vd_0 - vh_1))/2. + d*f_0*md_0*(-vd_0 + vh_1) + (d*f_1*mh_1*(-vd_0 + vh_1))/2. + (d*f_1*md_1*(-vd_1 + vh_2))/2. + d*(1 - f_0 - f_1)*mh_2*(-vd_1 + vh_2) + ((1 - d)*md_1*(-vh_1 + 2*vh_2))/2.))/2.));

double pdhtplus1=bound(pdh + 0.01*(f_0*(-((1 - d)*md_0*vd_0)/2. + ((1 - d)*rdd*md_0*(-vd_0 + vd_1))/2. + ((1 - d)*md_0*(-vd_0 + vh_1))/2. + d*f_0*md_0*(-vd_0 + vh_1) + (d*f_1*mh_1*(-vd_0 + vh_1))/2. + ((1 - d)*md_0*(-vd_0 + vd_1 + vh_1))/2. + (d*f_1*md_1*(-vd_1 + vh_2))/2. + d*(1 - f_0 - f_1)*mh_2*(-vd_1 + vh_2)) + (f_1*(-((1 - d)*mh_1*(2*vd_0 - vd_1))/2. + ((1 - d)*mh_1*vh_1)/2. + d*f_0*md_0*(-vd_0 + vh_1) + (d*f_1*mh_1*(-vd_0 + vh_1))/2. + ((1 - d)*md_1*(-vd_1 + vh_2))/2. + (d*f_1*md_1*(-vd_1 + vh_2))/2. + d*(1 - f_0 - f_1)*mh_2*(-vd_1 + vh_2) + ((1 - d)*rhd*md_1*(-vh_1 + vh_2))/2.))/2.));

gsl_vector_set(f,0, phhtplus1);

gsl_vector_set(f,1, pdhtplus1);

}


void write_params(void *params)
{
    // 
double d = ((struct rparams *) params)->d;
double md_0 = ((struct rparams *) params)->md_0;
double mh_1 = ((struct rparams *) params)->mh_1;
double md_1 = ((struct rparams *) params)->md_1;
double mh_2 = ((struct rparams *) params)->mh_2;
double pdhinit = ((struct rparams *) params)->pdhinit;
double phhinit = ((struct rparams *) params)->phhinit;
double f_0 = ((struct rparams *) params)->f_0;
double f_1 = ((struct rparams *) params)->f_1;
double pdh = ((struct rparams *) params)->pdh;
double phh = ((struct rparams *) params)->phh;
double rhh = ((struct rparams *) params)->rhh;
double rhd = ((struct rparams *) params)->rhd;
double rdd = ((struct rparams *) params)->rdd;
double vh_1 = ((struct rparams *) params)->vh_1;
double vh_2 = ((struct rparams *) params)->vh_2;
double vd_0 = ((struct rparams *) params)->vd_0;
double vd_1 = ((struct rparams *) params)->vd_1;



    // 
DataFile << endl << endl  << "d;" << d << endl
 << "md_0;" << md_0 << endl
 << "mh_1;" << mh_1 << endl
 << "md_1;" << md_1 << endl
 << "mh_2;" << mh_2 << endl
 << "pdhinit;" << pdhinit << endl
 << "phhinit;" << phhinit << endl
 << endl;
}


void write_data(void *params, int time)
{
    // 
double d = ((struct rparams *) params)->d;
double md_0 = ((struct rparams *) params)->md_0;
double mh_1 = ((struct rparams *) params)->mh_1;
double md_1 = ((struct rparams *) params)->md_1;
double mh_2 = ((struct rparams *) params)->mh_2;
double pdhinit = ((struct rparams *) params)->pdhinit;
double phhinit = ((struct rparams *) params)->phhinit;
double f_0 = ((struct rparams *) params)->f_0;
double f_1 = ((struct rparams *) params)->f_1;
double pdh = ((struct rparams *) params)->pdh;
double phh = ((struct rparams *) params)->phh;
double rhh = ((struct rparams *) params)->rhh;
double rhd = ((struct rparams *) params)->rhd;
double rdd = ((struct rparams *) params)->rdd;
double vh_1 = ((struct rparams *) params)->vh_1;
double vh_2 = ((struct rparams *) params)->vh_2;
double vd_0 = ((struct rparams *) params)->vd_0;
double vd_1 = ((struct rparams *) params)->vd_1;



    if (time < 0)
    {
        // 
DataFile << "time;vh_1;vd_1;vd_0;rdd;pdh;rhd;rhh;f_1;f_0;phh;vh_2;" << endl;
    }

    // 
DataFile << time << ";" << vh_1 << ";" << 
vd_1 << ";" << 
vd_0 << ";" << 
rdd << ";" << 
pdh << ";" << 
rhd << ";" << 
rhh << ";" << 
f_1 << ";" << 
f_0 << ";" << 
phh << ";" << 
vh_2 << ";" << 
 endl;
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
    gsl_vector *x_selgrad = gsl_vector_alloc(2);

    // initialize the struct with parameters
    struct rparams paramstruct; 
  
    // initialize command line argument things
    // see generate_cpp.py 
    // 
		paramstruct.d = atof(argv[2]);
		paramstruct.md_0 = atof(argv[3]);
		paramstruct.mh_1 = atof(argv[4]);
		paramstruct.md_1 = atof(argv[5]);
		paramstruct.mh_2 = atof(argv[6]);
		paramstruct.pdhinit = atof(argv[7]);
		paramstruct.phhinit = atof(argv[8]);
		paramstruct.f_0 = atof(argv[9]);
		paramstruct.f_1 = atof(argv[10]);
		paramstruct.pdh = atof(argv[11]);
		paramstruct.phh = atof(argv[12]);
		paramstruct.rhh = atof(argv[13]);
		paramstruct.rhd = atof(argv[14]);
		paramstruct.rdd = atof(argv[15]);
		paramstruct.vh_1 = atof(argv[16]);
		paramstruct.vh_2 = atof(argv[17]);
		paramstruct.vd_0 = atof(argv[18]);
		paramstruct.vd_1 = atof(argv[19]);

   

    // ranges for cycling: sometimes numerical iterations
    // won't resolve as solutions slightly cycle around
    // the final value. To see whether cyling behaviour
    // is repetitive over time, however, we create these
    // vectors that stores a series of past values of 
    // the values for the switching rates and compares
    // those to values that are currently found
    gsl_vector *phh_range = gsl_vector_alloc(10);
    gsl_vector *pdh_range = gsl_vector_alloc(10);

    for (int ik = 0; ik < 10; ++ik)
    {
        // initialize the vector
        gsl_vector_set(phh_range, ik, 0);        
        gsl_vector_set(pdh_range, ik, 0);
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

        bool condition_phh = (fabs(paramstruct.phh - gsl_vector_get(x_selgrad, 0)) < 1e-7 || paramstruct.phh >= 0.999) || paramstruct.phh <= 0.001;
        bool condition_pdh = (fabs(paramstruct.pdh - gsl_vector_get(x_selgrad, 1)) < 1e-7 || paramstruct.pdh >= 0.999) || paramstruct.pdh <= 0.001;

        if (condition_phh && condition_pdh)
        {
            paramstruct.phh = gsl_vector_get(x_selgrad, 0);
            paramstruct.pdh = gsl_vector_get(x_selgrad, 1);

            break;
        }
        paramstruct.phh = gsl_vector_get(x_selgrad, 0);
        paramstruct.pdh = gsl_vector_get(x_selgrad, 1);
        
        assert(isnan(paramstruct.phh) == 0);
        assert(isnan(paramstruct.pdh) == 0);
        
        if (iter > 50000)
        {
            bool found_in_range1 = false;
            bool found_in_range2 = false;

            bool done = false;

            for (int ik = 0; ik < 10; ++ik)
            {
                if (!found_in_range1 && fabs(gsl_vector_get(phh_range, ik)-paramstruct.phh) < 1e-10)
                {
                    found_in_range1 = true;
                }

                if (!found_in_range2 && fabs(gsl_vector_get(pdh_range, ik)-paramstruct.pdh) < 1e-10)
                {
                    found_in_range2 = true;
                }

                if (found_in_range1 && found_in_range2)
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
            gsl_vector_set(phh_range, ik, gsl_vector_get(phh_range, ik - 1));
            gsl_vector_set(pdh_range, ik, gsl_vector_get(pdh_range, ik - 1));
        }

        gsl_vector_set(phh_range, 0, paramstruct.phh);
        gsl_vector_set(pdh_range, 0, paramstruct.pdh);

        if (iter % 1 == 0)
        {
            write_data(&paramstruct,iter);
        }
    }

    write_data(&paramstruct,iter);
    write_params(&paramstruct);


    gsl_vector_free(x_p);
    gsl_vector_free(x_selgrad);
}

