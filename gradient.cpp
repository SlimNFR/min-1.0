//GRADIENT MODULE - SOURCE FILE//
////////////////////////////////

//C++ libraries
#include <vector>
#include <math.h>
#include <iostream>

//Code defined libraries
#include "gradient.h"
#include "globals.h"
#include "structure.h"
#include "utils.h"


//DEFINING THE MODULES' FUNCTIONS & INITIALISING THE EXTERNAL VARIABLES

namespace gradient {

//namespace variables

std::vector<double> exchange_x;
std::vector<double> exchange_y;
std::vector<double> exchange_z;

std::vector<double> anisotropy_x;
std::vector<double> anisotropy_y;
std::vector<double> anisotropy_z;



std::vector<double> zeeman_x;
std::vector<double> zeeman_y;
std::vector<double> zeeman_z;

std::vector<double> lagrangian_x;
std::vector<double> lagrangian_y;
std::vector<double> lagrangian_z;


std::vector<double> lambda_param; //array to store the gradient of the energy with respect to the lambda parameters (the constraints)

std::vector<double> total;
std::vector<double> real_energy;



//namespace functions

void total_f()
{

    total.resize(3*no_of_atoms+4,0); //give the "total" array a size

    real_energy.resize(3*no_of_atoms,0);


//    utils::compute_global_magnetisation();

    //call the gradient computing functions
    gradient::uniaxial_anisotropy_f(spin::x,spin::y,spin::z,spin::nx,spin::ny,spin::nz,spin::k,
                                     gradient::anisotropy_x,gradient::anisotropy_y,gradient::anisotropy_z); //gradient of the anis. term wrt to spins
    gradient::lambda_constraints_wrt_spins_f();//gradient of the constraints wrt each spin
    gradient::lambda_constraints_wrt_lambda_f();//gradient of the constraints wrt each constraint
    gradient::zeeman_f();


    switch (no_of_dimensions) {
    case 3:
        gradient::exchange3D_f(spin::x,spin::y,spin::z,gradient::exchange_x,gradient::exchange_y,gradient::exchange_z);
        break;
    case 2:
        gradient::exchange2D_f(spin::x,spin::y,spin::z,gradient::exchange_x,gradient::exchange_y,gradient::exchange_z);
        break;
    case 1:
        gradient::exchange1D_f(spin::x,spin::y,spin::z,gradient::exchange_x,gradient::exchange_y,gradient::exchange_z);
        break;
    default:
        std::cout<<"Ooops, something went wrong computing the total gradient";
        break;
    }


    //initialise the "total" array with the computed values

    for(int i = 0; i<no_of_atoms;i++)
        {
                gradient::total[i] = gradient::exchange_x[i]   +
                                     gradient::anisotropy_x[i] +
                                     gradient::zeeman_x[i]     +
                                     gradient::lagrangian_x[i];

                gradient::total[i+no_of_atoms] = gradient::exchange_y[i]   +
                                                 gradient::anisotropy_y[i] +
                                                 gradient::zeeman_y[i]     +
                                                 gradient::lagrangian_y[i];

                gradient::total[i+2*no_of_atoms] = gradient::exchange_z[i]   +
                                                   gradient::anisotropy_z[i] +
                                                   gradient::zeeman_z[i]     +
                                                   gradient::lagrangian_z[i];

                //std::cout<<i<<" "<<gradient::total[i]<<" "<<gradient::exchange_x[i]<<" "<<gradient::anisotropy_x[i]<<" "<<gradient::lagrangian_x[i]<<std::endl;
                //std::cout<<i<<" "<<gradient::total[i]<<" "<<gradient::exchange_y[i]<<" "<<gradient::anisotropy_y[i]<<" "<<gradient::lagrangian_y[i]<<std::endl;


                gradient::real_energy[i] = gradient::exchange_x[i]   +
                                           gradient::anisotropy_x[i] +
                                           gradient::zeeman_x[i];

                gradient::real_energy[i+no_of_atoms] = gradient::exchange_y[i]   +
                                                       gradient::anisotropy_y[i] +
                                                       gradient::zeeman_y[i];

                gradient::real_energy[i+2*no_of_atoms] = gradient::exchange_z[i]   +
                                                         gradient::anisotropy_z[i] +
                                                         gradient::zeeman_z[i];




        }


    for(int i = 0; i<4 ; i++)
    {
        gradient::total[i+3*no_of_atoms] = gradient::lambda_param[i];
    }



   // print the total grad.

//        for (unsigned int i = 0; i<gradient::total.size();i++)
//        {
//            std::cout<<gradient::total[i]<<" ";
//        }
//        std::cout<<std::endl;

}

void zeeman_f()
{

    gradient::zeeman_x.resize(no_of_atoms);
    gradient::zeeman_y.resize(no_of_atoms);
    gradient::zeeman_z.resize(no_of_atoms);


    for(int i = 0; i<no_of_atoms; i++)
    {
        gradient::zeeman_x[i] = - mu_s*H*hx;
        gradient::zeeman_y[i] = - mu_s*H*hy;
        gradient::zeeman_z[i] = - mu_s*H*hz;

    }





}


void lambda_constraints_wrt_lambda_f()
{

    gradient::lambda_param.resize(4);

    std::cout<<"Computing the gradient with respect to the lambda parameters"<<"\n";
    //gradient::lambda_param[0] = 0.0;
    gradient::lambda_param[0] = - no_of_atoms * spin::global[0]*( spin::global[2] - v0_constrain_z);
    gradient::lambda_param[1] = - no_of_atoms * spin::global[1]*( spin::global[2] - v0_constrain_z);
    gradient::lambda_param[2] = - no_of_atoms * spin::global[2]*( spin::global[2] - v0_constrain_z);

    gradient::lambda_param[3] = - (spin::global[2] - v0_constrain_z);//additional constraint to make the "v0 = - 0.5" case work


//    for (int i = 3; i<no_of_atoms+3;i++)
//    {

//    gradient::lambda_param[i] = - (sqrt(x[i-3]*x[i-3] + y[i-3]*y[i-3] + z[i-3]*z[i-3]) - 1.0);


//    }



}


void lambda_constraints_wrt_spins_f()
{

    //function that computes the gradient of the orientation constraint term with respect to each spin

    std::cout<<"Computing the gradient of the orientation constraint term"<<"\n";



    gradient::lagrangian_x.resize(no_of_atoms,0);
    gradient::lagrangian_y.resize(no_of_atoms,0);
    gradient::lagrangian_z.resize(no_of_atoms,0);



    for(int i = 0; i < no_of_atoms; i++)
    {

        gradient::lagrangian_x[i] = no_of_atoms*spin::lambda[0]*(spin::global[2]-v0_constrain_z);


        gradient::lagrangian_y[i] = no_of_atoms*spin::lambda[1]*(spin::global[2]-v0_constrain_z);

        gradient::lagrangian_z[i] = no_of_atoms*spin::lambda[2]*(2*spin::global[2] - v0_constrain_z) + spin::lambda[3];


    }



}


void uniaxial_anisotropy_f(std::vector<double> x, std::vector<double> y, std::vector<double> z, std::vector<double> nx, std::vector<double> ny, std::vector<double> nz,
                           std::vector<double>& k, std::vector<double>& g_anis_x, std::vector<double>& g_anis_y, std::vector<double>& g_anis_z)
{

   g_anis_x.resize(no_of_atoms,0);
   g_anis_y.resize(no_of_atoms,0);
   g_anis_z.resize(no_of_atoms,0);

   switch (uniform_init_anisotropy) {
   case 1:
    std::cout<<"Computing the anis. gradient: uniform case."<<std::endl;
    for(int i = 0;i < no_of_atoms;i++)
    {

        g_anis_x[i] = - 2 * k[0] * (x[i]*nx[i] + y[i]*ny[i] + z[i]*nz[i]) * nx[i];
        g_anis_y[i] = - 2 * k[0] * (x[i]*nx[i] + y[i]*ny[i] + z[i]*nz[i]) * ny[i];
        g_anis_z[i] = - 2 * k[0] * (x[i]*nx[i] + y[i]*ny[i] + z[i]*nz[i]) * nz[i];


    }



       break;
   case 0:

       std::cout<<"Computing the anis. gradient: non-uniform case."<<std::endl;
       for(int i = 0; i < no_of_atoms/2; i++)
       {
           g_anis_x[i] = - 2 * k[0] * (x[i]*nx[i] + y[i]*ny[i] + z[i]*nz[i]) * nx[i];
           g_anis_y[i] = - 2 * k[0] * (x[i]*nx[i] + y[i]*ny[i] + z[i]*nz[i]) * ny[i];
           g_anis_z[i] = - 2 * k[0] * (x[i]*nx[i] + y[i]*ny[i] + z[i]*nz[i]) * nz[i];

       }

       for(int i = no_of_atoms/2; i < no_of_atoms; i++)
       {

           g_anis_x[i] = - 2 * k[1] * (x[i]*nx[i] + y[i]*ny[i] + z[i]*nz[i]) * nx[i];
           g_anis_y[i] = - 2 * k[1] * (x[i]*nx[i] + y[i]*ny[i] + z[i]*nz[i]) * ny[i];
           g_anis_z[i] = - 2 * k[1] * (x[i]*nx[i] + y[i]*ny[i] + z[i]*nz[i]) * nz[i];

       }

       break;
   default:
       break;
   }


}




void exchange3D_f(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z,
                  std::vector<double>& g_exch_x, std::vector<double>& g_exch_y, std::vector<double>& g_exch_z)
{
    try
    {

        if(no_of_dimensions != 3)throw no_of_dimensions;


    }
    catch(const int x){std::cout<<"no. of dimensions is incorrrectly set, therefore the 3D exchange gradient cannot be computed"<<"\n";
                       return;}

    std::cout<<"Computing the exchange gradient: 3D case"<<"\n";

    g_exch_x.resize(no_of_atoms,0);
    g_exch_y.resize(no_of_atoms,0);
    g_exch_z.resize(no_of_atoms,0);


    switch (uniform_exchange_interaction) {
    case 1:
    {
        for(int k = 0; k < no_of_atoms; k = k + n1*n2)//the horizontal interactions
        {
        for(int i = 0; i < n1*n2;i = i + n2)
        {
            g_exch_x[i + k] = - J * x[i + k + 1]; //
            g_exch_y[i + k] = - J * y[i + k + 1]; // initialise the first three values of the g_exch arrays with the components of the second spin on each line
            g_exch_z[i + k] = - J * z[i + k + 1]; //

            g_exch_x[i + k + n2 - 1] = - J * x[i + k + n2 - 2];
            g_exch_y[i + k + n2 - 1] = - J * y[i + k + n2 - 2]; //initialise the last three values of the g_exch array with the components of the (N-1)th spin on each line
            g_exch_z[i + k + n2 - 1] = - J * z[i + k + n2 - 2];

            for(int j = 1;j <= n2 - 2;j++)
            {
                g_exch_x[i+j+k] = - J * (x[i+j+k-1] + x[i+j+k+1]);
                g_exch_y[i+j+k] = - J * (y[i+j+k-1] + y[i+j+k+1]);
                g_exch_z[i+j+k] = - J * (z[i+j+k-1] + z[i+j+k+1]);


            }
        }
        }

        for(int k = 0; k < no_of_atoms; k = k + n1*n2)
        {
        for(int i = 0; i < n2;i = i + 1)//the vertical interactions
        {
            g_exch_x[i + k] = g_exch_x[i + k] - J * x[i + k + n2]; //
            g_exch_y[i + k] = g_exch_y[i + k] - J * y[i + k + n2]; // initialise the first three values of the g_exch arrays with the components of the second spin on each line
            g_exch_z[i + k] = g_exch_z[i + k] - J * z[i + k + n2]; //







            g_exch_x[i + k + n1*n2 - n2] = g_exch_x[i + k + n1*n2 - n2] - J * (x[i + k + n1*n2 - 2*n2]);
            g_exch_y[i + k + n1*n2 - n2] = g_exch_y[i + k + n1*n2 - n2] - J * (y[i + k + n1*n2 - 2*n2]); //initialise the last three values of the g_exch array with the components of the (N-1)th spin on each line
            g_exch_z[i + k + n1*n2 - n2] = g_exch_z[i + k + n1*n2 - n2] - J * (z[i + k + n1*n2 - 2*n2]);


            for(int j = n2;j < n1*n2 - n2;j = j + n2)
            {
                g_exch_x[i+j+k] = g_exch_x[i+j+k] - J * ( x[i+j+k-n2] + x[i+j+k+n2] );
                g_exch_y[i+j+k] = g_exch_y[i+j+k] - J * ( y[i+j+k-n2] + y[i+j+k+n2] );
                g_exch_z[i+j+k] = g_exch_z[i+j+k] - J * ( z[i+j+k-n2] + z[i+j+k+n2] );


            }
        }
        }


                for(int i = 0; i < n2; i = i + 1) //the in-depth interactions
                {
                for(int j = 0; j < n1*n2;j = j + n2)
                {
                    g_exch_x[i + j] = g_exch_x[i + j] - J * x[i + j + n1*n2]; //
                    g_exch_y[i + j] = g_exch_y[i + j] - J * y[i + j + n1*n2]; // initialise the first three values of the g_exch arrays with the components of the second spin on each line
                    g_exch_z[i + j] = g_exch_z[i + j] - J * z[i + j + n1*n2];







                    g_exch_x[i + j + no_of_atoms - n1*n2] = g_exch_x[i + j + no_of_atoms - n1*n2] - J * (x[i + j + no_of_atoms - 2*n1*n2]);
                    g_exch_y[i + j + no_of_atoms - n1*n2] = g_exch_y[i + j + no_of_atoms - n1*n2] - J * (y[i + j + no_of_atoms - 2*n1*n2]); //initialise the last three values of the g_exch array with the components of the (N-1)th spin on each line
                    g_exch_z[i + j + no_of_atoms - n1*n2] = g_exch_z[i + j + no_of_atoms - n1*n2] - J * (z[i + j + no_of_atoms - 2*n1*n2]);


                    for(int k = n1*n2;k < no_of_atoms - n1*n2;k = k + n1*n2)
                    {
                        g_exch_x[i+j+k] = g_exch_x[i+j+k] - J * ( x[i+j+k-n1*n2] + x[i+j+k+n1*n2] );
                        g_exch_y[i+j+k] = g_exch_y[i+j+k] - J * ( y[i+j+k-n1*n2] + y[i+j+k+n1*n2] );
                        g_exch_z[i+j+k] = g_exch_z[i+j+k] - J * ( z[i+j+k-n1*n2] + z[i+j+k+n1*n2] );


                    }
                }
                }




        break;
    }
    case 0:
    {


    break;
    }

    default:
        break;
    }



}




void exchange2D_f(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z,
                  std::vector<double>& g_exch_x, std::vector<double>& g_exch_y, std::vector<double>& g_exch_z)
{
    try
    {

        if(no_of_dimensions != 2)throw no_of_dimensions;


    }
    catch(const int x){std::cout<<"no.of dimensions is incorrectly set, therefore the 2D exchange gradient cannot be computed"<<"\n";
                       return;}

    std::cout<<"Computing the exchange gradient: 2D case"<<"\n";

    g_exch_x.resize(no_of_atoms,0);
    g_exch_y.resize(no_of_atoms,0);
    g_exch_z.resize(no_of_atoms,0);


    switch (uniform_exchange_interaction) {
    case 1:
    {







    for(int i = 0; i < no_of_atoms;i = i + n2)//the horizontal interactions
    {
        g_exch_x[i] = - J * x[i+1]; //
        g_exch_y[i] = - J * y[i+1]; // initialise the first three values of the g_exch arrays with the components of the second spin on each line
        g_exch_z[i] = - J * z[i+1]; //

        g_exch_x[i+n2 - 1] = - J * x[i+n2 - 2];
        g_exch_y[i+n2 - 1] = - J * y[i+n2 - 2]; //initialise the last three values of the g_exch array with the components of the (N-1)th spin on each line
        g_exch_z[i+n2 - 1] = - J * z[i+n2 - 2];

        for(int j = 1;j <= n2 - 2;j++)
        {
            g_exch_x[i+j] = - J * (x[i+j-1] + x[i+j+1]);
            g_exch_y[i+j] = - J * (y[i+j-1] + y[i+j+1]);
            g_exch_z[i+j] = - J * (z[i+j-1] + z[i+j+1]);


        }
    }


    for(int i = 0; i < n2;i = i + 1)//the vertical interactions
    {
        g_exch_x[i] = g_exch_x[i] - J * x[i + n2]; //
        g_exch_y[i] = g_exch_y[i] - J * y[i + n2]; // initialise the first three values of the g_exch arrays with the components of the second spin on each line
        g_exch_z[i] = g_exch_z[i] - J * z[i + n2]; //







        g_exch_x[i + no_of_atoms - n2] = g_exch_x[i + no_of_atoms - n2] - J * (x[i + no_of_atoms - 2*n2]);
        g_exch_y[i + no_of_atoms - n2] = g_exch_y[i + no_of_atoms - n2] - J * (y[i + no_of_atoms - 2*n2]); //initialise the last three values of the g_exch array with the components of the (N-1)th spin on each line
        g_exch_z[i + no_of_atoms - n2] = g_exch_z[i + no_of_atoms - n2] - J * (z[i + no_of_atoms - 2*n2]);


        for(int j = n2;j < no_of_atoms - n2;j = j + n2)
        {
            g_exch_x[i+j] = g_exch_x[i+j] - J * ( x[i+j-n2] + x[i+j+n2] );
            g_exch_y[i+j] = g_exch_y[i+j] - J * ( y[i+j-n2] + y[i+j+n2] );
            g_exch_z[i+j] = g_exch_z[i+j] - J * ( z[i+j-n2] + z[i+j+n2] );


        }
    }
    break;
    }


    case 0:
        std::cout<<"2D: Non-uniform exchange interaction not yet implemented  --- couldn't compute its gradient"<<std::endl;
        break;

    default:
        std::cout<<"Ooops, something went wrong trying to compute the gradient of the ex. energy (1D)"<<std::endl;
        break;

}
}





void exchange1D_f(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z,
                  std::vector<double>& g_exch_x, std::vector<double>& g_exch_y, std::vector<double>& g_exch_z)
{

    try
    {

        if(no_of_dimensions == 3)throw no_of_dimensions;


    }
    catch(const int x){std::cout<<"no.of dimensions is incorrectly set, therefore the 1D exchange gradient cannot be computed"<<"\n";
                       return;}


    g_exch_x.resize(no_of_atoms);
    g_exch_y.resize(no_of_atoms);
    g_exch_z.resize(no_of_atoms);

    std::cout<<"Computing the exchange gradient: 1D case"<<"\n";

    switch (uniform_exchange_interaction){

    case 1:
    {

        std::cout<<"bla";

        g_exch_x[0] = - J * x[1]; //
        g_exch_y[0] = - J * y[1]; // initialise the first three values of the g_exch arrays with the components of the second spin
        g_exch_z[0] = - J * z[1]; //

        g_exch_x[no_of_atoms - 1] = - J * x[no_of_atoms - 2];
        g_exch_y[no_of_atoms - 1] = - J * y[no_of_atoms - 2]; //initialise the last three values of the g_exch array with the components of the (N-1)th spin
        g_exch_z[no_of_atoms - 1] = - J * z[no_of_atoms - 2];


        for(int i = 1; i <= no_of_atoms - 2; i++)
        {

            g_exch_x[i] = - J * (x[i-1] + x[i + 1]);
            g_exch_y[i] = - J * (y[i-1] + y[i + 1]); //fill the rest of the array; besides the first and the last but one spin, all the others occur two times in the derivative of the exchange
            g_exch_z[i] = - J * (z[i-1] + z[i + 1]);

        }
        std::cout<<"end";

        break;

    }
    case 0:
    {
        std::cout<<"1D: Non-uniform exchange interaction not yet implemented  --- couldn't compute its gradient"<<std::endl;
        break;
    }



    default:
        std::cout<<"Ooops, something went wrong trying to compute the gradient of the ex. energy (1D)"<<std::endl;
        break;


    }
}


}
