//ENERGY MODULE - SOURCE FILE//
//////////////////////////////


//C++ libraries
#include <vector>
#include <math.h>
#include <iostream>

//Code defined libraries
#include "globals.h"
#include "energy.h"
#include "utils.h"
#include "structure.h"



//DEFINING THE MODULES' FUNCTIONS & INITIALISING THE EXTERNAL VARIABLES

namespace energy{

//namespace variables
double exchange;
double uniaxial_anisotropy;
double lagrangian_term;
double total;
double real_energy;
double zeeman;


//namespace functions

void total_f()
{

    energy::uniaxial_anisotropy = energy::uniaxial_anisotropy_f(spin::nx,spin::ny,spin::nz,
                                                                spin::x,spin::y,spin::z,
                                                                spin::k);
    energy::lagrangian_term = energy::lagrangian_term_f(spin::x,spin::y,spin::z,
                                                        spin::lambda);

    energy::exchange = energy::exchange_f(spin::x,spin::y,spin::z);

    energy::zeeman = energy::zeeman_f();

    energy::total = energy::uniaxial_anisotropy + energy::exchange + energy::zeeman + energy::lagrangian_term;

    energy::real_energy = energy::uniaxial_anisotropy + energy::exchange + energy::zeeman;






    std::cout<<"Anisotropy energy: "<<energy::uniaxial_anisotropy<<" J"<<std::endl;
    std::cout<<"Exchange energy: " <<energy::exchange<<" J"<<std::endl;
    std::cout<<"Lagrang energy: " <<energy::lagrangian_term<<" J"<<std::endl;
    std::cout<<"Zeeman energy: " <<energy::zeeman<<" J"<<std::endl;
    std::cout<<"TOTAL ENERGY OF THE SYSTEM IS: "<<energy::total<<" J"<<"\n\n";
    std::cout<<"Real energy of the system is: "<<energy::real_energy<<" J"<<"\n\n";


}

double zeeman_f()
{

    energy::zeeman = 0.0;

    for(int i = 0; i<no_of_atoms; i++)
    {
        energy::zeeman = energy::zeeman - (spin::x[i]*hx +
                                           spin::y[i]*hy +
                                           spin::z[i]*hz ) * mu_s*H;

    }

    return energy::zeeman;

}

double lagrangian_term_f(std::vector<double> x, std::vector<double> y, std::vector<double> z,
                         std::vector<double> lambda)
{


    double lagrang_term = 0; // in this variable, we will store the sum of the lagrangian terms

    utils::compute_global_magnetisation(); // call to function which computes the global magn. vector

    lagrang_term = lagrang_term + no_of_atoms * ( lambda[0] * ( spin::global[0] - v0_constrain_x) +
                                                  lambda[1] * ( spin::global[1] - v0_constrain_y) +
                                                  lambda[2] * ( spin::global[2] - v0_constrain_z) ); //compute the vectorial constraint
    //    for(int i = 0;i < no_of_atoms ;i++)
    //    {


    //        lagrang_term = lagrang_term + lambda[i+3] * ( sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]) - 1.0); //compute the modulus constraints
    //                                                                                                                         //the i+3 index refers to the fact that the first 3
    //                                                                                                                         // components of the lambda array belong to the first constraint, which is a vectorial quantity

    //    }

    return lagrang_term;
}


double exchange_f(std::vector<double> x, std::vector<double> y, std::vector<double> z)
{


    double ex_energy = 0;



    switch(no_of_dimensions)
    {

    case 1:
    {
        std::cout<<"1-dimensional case: computing the exchange energy .."<<std::endl;
        for(int i = 0; i<no_of_atoms-1; i++)
        {
            ex_energy = ex_energy - J * ( x[i]*x[i+1] +
                    y[i]*y[i+1] +
                    z[i]*z[i+1] ); //compute the i--->i+1 interaction


        }

        break;
    }
    case 2:
    {
        std::cout<<"2-dimensional case: computing the exchange energy .."<<std::endl;



        for(int i = 0; i < no_of_atoms;i = i + n2)
        {

            for(int j = 0;j < n2 - 1;j++)
            {
                ex_energy = ex_energy - J * ( x[i+j]*x[i+j+1] +
                        y[i+j]*y[i+j+1] +
                        z[i+j]*z[i+j+1]);//compute the horizontal interactions

            }
        }


        for(int i = 0; i < no_of_atoms - n2; i = i + n2)
        {

            for(int j = 0; j < n2; j++)
            {

                ex_energy = ex_energy - J * ( x[i+j]*x[i+j+n2] +
                        y[i+j]*y[i+j+n2] +
                        z[i+j]*z[i+j+n2] );//compute the vertical interactions

            }

        }




        break;
    }
    case 3:
    {    std::cout<<"3-dimensional case: computing the exchange energy .."<<std::endl;

        //compute the horizontal interactions
        for(int k = 0; k < no_of_atoms; k = k + n2*n1)
        {
            for(int i = 0; i< n2*n1; i = i + n2 )
            {
                for(int j = 0; j < n2 - 1; j = j + 1)
                {
                    ex_energy = ex_energy - J * ( x[i+j+k]*x[i+j+k+1] +
                            y[i+j+k]*y[i+j+k+1] +
                            z[i+j+k]*z[i+j+k+1] );
                    std::cout<<ex_energy<<std::endl;
                }
            }
        }

        //compute the vertical interactions
        for(int k = 0; k < no_of_atoms; k = k + n2*n1)
        {
            for(int i = 0; i< n2*n1 - n2 ; i = i + n2 )
            {
                for(int j = 0; j < n2 ; j = j + 1)
                {
                    ex_energy = ex_energy - J * ( x[i+j+k]*x[i+j+k+n2] +
                            y[i+j+k]*y[i+j+k+n2] +
                            z[i+j+k]*z[i+j+k+n2] );
                    std::cout<<ex_energy<<std::endl;
                }
            }
        }


        //compute the in-depth interactions
        for(int k = 0; k < no_of_atoms - n2*n1; k = k + n2*n1)
        {
            for(int i = 0; i< n2*n1; i = i + n2 )
            {
                for(int j = 0; j < n2 ; j = j + 1)
                {
                    ex_energy = ex_energy - J * ( x[i+j+k]*x[i+j+k+n1*n2] +
                            y[i+j+k]*y[i+j+k+n1*n2] +
                            z[i+j+k]*z[i+j+k+n1*n2] );

                    std::cout<<ex_energy<<std::endl;
                }
            }
        }

        break;
    }
    default:
        std::cout<<"Ooops, couldn't compute the exchange energy."<<std::endl;


    }







    return ex_energy;

}


double uniaxial_anisotropy_f(std::vector<double> nx, std::vector<double> ny, std::vector<double> nz,
                             std::vector<double> mx, std::vector<double> my, std::vector<double> mz,
                             std::vector<double> k)
{
    double anis_energy = 0;


    switch (uniform_init_anisotropy) {
    case 1: //all the spins have the same anisotropy constant and the same orientation of the uni-axial vector
    {

        for(int i = 0; i<no_of_atoms; i++)
        {


            anis_energy = anis_energy - k[0] * ( pow( mx[i]*nx[i] + my[i]*ny[i] + mz[i]*nz[i], 2.0));


        }



        break;
    }
    case 0://half the spins have the same anisotropy constant and the same orientation of the uni-axial vector
    {
        for(int i = 0; i<no_of_atoms/2; i++)
        {


            anis_energy = anis_energy - k[0] * ( pow( mx[i]*nx[i] + my[i]*ny[i] + mz[i]*nz[i], 2));




        }


        for(int i = no_of_atoms/2; i<no_of_atoms; i++)
        {


            anis_energy = anis_energy - k[1] * ( pow( mx[i]*nx[i] + my[i]*ny[i] + mz[i]*nz[i], 2));




        }
        break;

    }
    default:
        std::cout<<"Oooppss, something went wrong computing the uni-axial anisotropy energy"<<std::endl;
        break;
    }





    return anis_energy;

}
}
