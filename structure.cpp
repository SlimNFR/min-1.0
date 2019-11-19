//STRUCTURE MODULE - SOURCE FILE//
//////////////////////////////


//C++ libraries
#include <vector>
#include <math.h>
#include <iostream>

//Code defined libraries
#include "structure.h"
#include "globals.h"
#include "energy.h"
#include "utils.h"


//DEFINING THE MODULES' FUNCTIONS & INITIALISING THE EXTERNAL VARIABLES

namespace spin {

//namespace variables

std::vector<double> x; //saves the x component of spins
std::vector<double> y; //saves the y component of spins
std::vector<double> z; //saves the z component of spins


std::vector<double> nx; //saves the x component of the uni-axial vector for all the spins
std::vector<double> ny; //saves the y component of the uni-axial vector for all the spins
std::vector<double> nz; //saves the z component of the uni-axial vector for all the spins

std::vector<double> k;  //saves the anisotropy constant for each spin


std::vector<double> lambda; //saves the lagrangian constraints

std::vector<double> global={0,0,0};

//namespace functions


void print_f()
{
    std::cout<<"Printing the spins' components:"<<"\n";
    for(int i = 0; i<no_of_atoms;i++)
    {
        std::cout<<spin::x[i]<<","<<spin::y[i]<<","<<spin::z[i]<<"\n";
    }

}

}

namespace structure {


//namespace variables



//namespace functions

void initial_structure_f(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z)
{

    std::cout<<"GENERATING INITIAL "<<n1<<"x"<<n2<<"x"<<n3<<" SPIN STRUCTURE "<<std::endl;
    std::cout<<"TOTAL NO. OF SPINS: "<<no_of_atoms<<std::endl<<std::endl;




    switch (uniform_init_structure) {
    case 1:
    {
        std::cout<<"UNIFORM SPIN ORIENTATION"<<std::endl<<std::endl;
        x.resize(no_of_atoms,1);
        y.resize(no_of_atoms,0);
        z.resize(no_of_atoms,0);
        break;
    }

    case 0:
    {
        std::cout<<"NON-UNIFORM SPIN ORIENTATION"<<std::endl<<std::endl;
        x.resize(no_of_atoms);
        y.resize(no_of_atoms);
        z.resize(no_of_atoms);
        for(int i = 0; i < no_of_atoms/2; i++)
        {
            x[i] = 1;
            y[i] = 0;      //spins on z-axis
            z[i] = 0;
        }

        for(int i = no_of_atoms/2; i < no_of_atoms; i++)
        {
            x[i] = 0;
            y[i] = 1;      //spins on y-axis
            z[i] = 0;
        }
    break;
    }
    case -1:
    {
        std::cout<<"GIVE ME THE SPIN COMPONENTS"<<"\n";

        x.resize(no_of_atoms);
        y.resize(no_of_atoms);
        z.resize(no_of_atoms);

        for(int i = 0;i<no_of_atoms;i++)
        {
             std::cout<<"x["<<i<<"],"<<"y["<<i<<"],"<<"z["<<i<<"]"<<"\n";
             std::cin>>x[i];std::cin>>y[i];std::cin>>z[i];

        }


    break;
    }

    default:
        break;
    }





}


void uniaxial_anisotropy_f(std::vector<double>& nx, std::vector<double>& ny, std::vector<double>& nz,
                           std::vector<double>& k)
{



    switch (uniform_init_anisotropy) {
    case 1: //if the uniform flag is set, the anisotropy vector has the same
    {        //orientation for each spin; the anis. constant is also identical.
        std::cout<<"Setting up: homogenous uniaxial anisotropy"<<std::endl;

        nx.resize(no_of_atoms,0);
        ny.resize(no_of_atoms,1.0);
        nz.resize(no_of_atoms,0);
        k.resize(1,0.38*1e-1); // an array with only one element : the unique anisotropy constant



        //Fe: 0.801*1e-4
        //FePt: 0.38*1e-1
        break;
    }
    case 0:  //the anisotropy vector and the magnitude vary from spin to spin
    {         // TEMPORARY: half the spins have k1, the others k2

        std::cout<<"Setting up: inhomogenous uniaxial anisotropy"<<std::endl;

        nx.resize(no_of_atoms,0);
        ny.resize(no_of_atoms,0);
        nz.resize(no_of_atoms,0);
        k.resize(2,0); // let's assume two different anisotropy constants for now


        for(int i = 0; i < no_of_atoms/2; i++)
        {
            nx[i] = 0.0;
            ny[i] = 0.0;      //anisotropy on z-axis
            nz[i] = 1.0;
        }

        for(int i = no_of_atoms/2; i < no_of_atoms; i++)
        {
            nx[i] = 0.0;
            ny[i] = 1.0;      //anisotropy on y-axis
            nz[i] = 0.0;
        }

        k[0] = 0.801*1e-4;   //anisotropy
        k[1] = 0.801*1e-4;   //constants


        break;
    }
    default:
        std::cout<<"Ooops, something went wrong setting setting up the uniaxial anisotropy"<<std::endl;
        break;
    }


}



void init_lagrangian_constraints_f(std::vector<double>& lambda)
{

    lambda.resize(4);

    switch (uniform_lagrangian_constraints) {
    case 1:
    {
        std::cout<<"Setting up the uniform lagrangian constraints"<<std::endl;
        lambda[0] = lambda_x;
        lambda[1] = lambda_y;  //setup the vectorial components of the lagrangian constraints
        lambda[2] = lambda_z;

        lambda[3] = l1;

//        for(int i = 3;i < no_of_atoms + 3; i++) //setup the rest of the constraints
//        {

//            lambda[i] = 1;

//        }

        break;
    }

    case -1://user input
    {
        std::cout<<"GIVE ME THE LAMBDA PARAMS:"<<"\n";

        for(int i = 0; i< 3;i++)
            std::cin>>lambda[i];


        break;
    }
    default:
        std::cout<<"Ooops, something went wrong setting up the lagrangian constraints"<<std::endl;
        break;
    }


}


}
