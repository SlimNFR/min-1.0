//GRADIENT MODULE - SOURCE FILE//
////////////////////////////////

#ifndef GRADIENT_H
#define GRADIENT_H

//C++ libraries
#include <vector>
#include <iostream>

//Code defined libraries
#include "globals.h"


//namespace for the gradient terms
namespace gradient {

//namespace variables

extern std::vector<double> exchange_x; //
extern std::vector<double> exchange_y; // arrays to store the gradient of the exchange energy term with respect to the spins' components on x,y,z
extern std::vector<double> exchange_z; //

extern std::vector<double> anisotropy_x;
extern std::vector<double> anisotropy_y; //arrays to store the gradient of the anis.energy term with respect to the spins' components on x, y, z
extern std::vector<double> anisotropy_z;

extern std::vector<double> zeeman_x;
extern std::vector<double> zeeman_y; //arrays to store the gradient of the zeeman energy with respect to the spins' compontents on x,y,z
extern std::vector<double> zeeman_z;

extern std::vector<double> lagrangian_x;
extern std::vector<double> lagrangian_y; //arrays to store the gradient of the lagrangian term with respect to the spins' components on x,y,z
extern std::vector<double> lagrangian_z;


extern std::vector<double> lambda_param; //array to store the gradient of the energy with respect to the lambda parameters (the constraints)


extern std::vector<double> total; //array to store the total gradient
extern std::vector<double> real_energy; //array to store the gradient of the real energy

//namespace functions

void total_f();

void zeeman_f();

void lambda_constraints_wrt_lambda_f();//function which computes the gradient of the lagrangian terms with respect to the lambda parameters



void lambda_constraints_wrt_spins_f(); //function which computes the gradient of the lagrangian terms with respect to the spins' components

void uniaxial_anisotropy_f(std::vector<double> x, std::vector<double> y, std::vector<double> z,
                           std::vector<double> nx, std::vector<double> ny, std::vector<double> nz,
                           std::vector<double>& k,
                           std::vector<double>& g_anis_x, std::vector<double>& g_anis_y, std::vector<double>& g_anis_z); //function which computes the gradient of the uniaxial anisotropy term


void exchange3D_f(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z,
                  std::vector<double>& g_exch_x, std::vector<double>& g_exch_y, std::vector<double>& g_exch_z); //computes the gradient of the exchange energy in the 3D case


void exchange2D_f(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z,
                  std::vector<double>& g_exch_x, std::vector<double>& g_exch_y, std::vector<double>& g_exch_z); //computes the gradient of the exchange energy in the 2D case

void exchange1D_f(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z,
                  std::vector<double>& g_exch_x, std::vector<double>& g_exch_y, std::vector<double>& g_exch_z); //computes the gradient of the exchange energy in the 1D case

}




#endif // GRADIENT_H
