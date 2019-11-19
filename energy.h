//ENERGY MODULE - HEADER FILE//
//////////////////////////////



#ifndef ENERGY_H
#define ENERGY_H



//C++ libraries
#include<vector>


//Code defined libraries


//namespace for the energy terms which need to be minimised
namespace energy {


//namespace variables

extern double exchange; //variable to save the exchange energy
extern double uniaxial_anisotropy; //variable to save the anisotropy energy
extern double lagrangian_term; //variable to save the energy corresponding to the lagrangian terms
extern double total; //variable to save the total energy of the system
extern double real_energy; //variable to save the real energy of the system
extern double zeeman;


//namespace functions

void total_f(); //func to compute the total energy of the system

double zeeman_f();//func to compute the zeeman energy

double lagrangian_term_f(std::vector<double> x, std::vector<double> y, std::vector<double> z,
                         std::vector<double> lambda); //function to compute the energy related to the lagrangian terms

double exchange_f(std::vector<double> x, std::vector<double> y, std::vector<double> z);//function to compute the exchange energy: 1D, 2D & 3D case

double uniaxial_anisotropy_f(std::vector<double> nx, std::vector<double> ny, std::vector<double> nz,
                             std::vector<double> mx, std::vector<double> my, std::vector<double> mz,
                             std::vector<double> k); //function to calculate the uniaxial anisotropy energy






}


#endif // ENERGY_H
