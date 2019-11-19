 //UTILS MODULE - SOURCE FILE//
//////////////////////////////


//C++ libraries
#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

//Code defined libraries
#include "utils.h"
#include "globals.h"
#include "structure.h"
#include "gradient.h"



//DEFINING THE MODULES' FUNCTIONS & INITIALISING THE EXTERNAL VARIABLES


namespace utils {


//namespace variables

std::vector<double>global_magnetisation = {0,0,0};
std::vector<double> sums = {0,0,0,0};



std::vector<double>temp_x;
std::vector<double>temp_y;
std::vector<double>temp_z;
std::vector<double>temp_lambda;
std::vector<double>temp_gradient;
double temp_energy;




//namespace functions
double compute_spin_magnitude(double x, double y, double z)
{


return sqrt(x*x + y*y + z*z);


}



void print_all(std::ofstream& f1, std::ofstream& f2,int step)
{




for(int i = 0;i<no_of_atoms;i++)
{


   f1<<std::setprecision(10);
   f1<<i<<" "<<spin::x[i]<<" "<<spin::y[i]<<" "<<spin::z[i]<<" "<<sqrt(spin::x[i]*spin::x[i]+spin::y[i]*spin::y[i]+spin::z[i]*spin::z[i])<<" "<<step<<"\n";


//    <<spin::lambda[i]<<" "<<spin::lambda[i+1]<<" "<<spin::lambda[i+2]<<" "<<" "<<spin::lambda[i+3] << "  "
//   <<gradient::total[i]<<" "<<gradient::total[i+no_of_atoms]<<" "<<gradient::total[i+2*no_of_atoms]<<" "<<gradient::total[i+3*no_of_atoms + 3]<<"\n";

}




f2<<std::setprecision(10);
f2<<energy::real_energy<<" "<<utils::compute_modul(gradient::real_energy)<<" "<<step<<"\n";


    /*<<gradient::total[*no_of_atoms]<<" "<<gradient::total[4*no_of_atoms + 1]<<gradient::total[4*no_of_atoms + 2];*/



}



void save_state(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z,
                std::vector<double>& lambda,
                std::vector<double>& gradient,
                double& energy)
{

x.resize(no_of_atoms);
y.resize(no_of_atoms);
z.resize(no_of_atoms);
lambda.resize(no_of_atoms+3);
gradient.resize(4*no_of_atoms+3);


x = spin::x;
y = spin::y;
z = spin::z;
lambda = spin::lambda;
gradient = gradient::total;
energy = energy::total;




}


double compute_modul(std::vector<double> vec)
{
    double sum = 0;

    for(unsigned int i = 0;i<vec.size();i++)
    {
        sum = sum + vec[i]*vec[i];

    }

    return sqrt(sum);
}


void compute_global_magnetisation()
{
    spin::global.resize(3);

    double modulus = 0, term1 = 0, term2 = 0, term3 = 0;

    for(int i = 0;i < no_of_atoms;i++)
    {


    term1 = term1 + spin::x[i];
    term2 = term2 + spin::y[i];
    term3 = term3 + spin::z[i];

    }

    modulus = sqrt( pow(term1,2.0) + pow(term2,2.0) + pow(term3,2.0)  );

    spin::global[0] = term1/modulus;
    spin::global[1] = term2/modulus;
    spin::global[2] = term3/modulus;


    std::cout<<"global magnetisation is: ("<<spin::global[0]<<","<<spin::global[1]<<","<<spin::global[2]<<")"<<std::endl; // ---- uncomment this if you want to print the global magn


}

void compute_magnetisation_sums(std::vector<double> x, std::vector<double> y, std::vector<double> z,
                                std::vector<double>& sums)
{


    double modulus = 0, term1 = 0, term2 = 0, term3 = 0;

    for(int i = 0;i < no_of_atoms;i++)
    {

    term1 = term1 + x[i];
    term2 = term2 + y[i];
    term3 = term3 + z[i];

    }

    modulus =  pow(term1,2.0) + pow(term2,2.0) + pow(term3,2.0); // sum of the x components squared + sum of the y components squared + sum of the z components squared

    sums[0] = term1;
    sums[1] = term2;
    sums[2] = term3;
    sums[3] = modulus;



}




}
