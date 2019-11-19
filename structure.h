//STRUCTURE MODULE - SOURCE FILE//
//////////////////////////////


#ifndef STRUCTURE_H
#define STRUCTURE_H


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


//namespace for spins related variables and functions
namespace spin {

//namespace variables

extern std::vector<double> x; //saves the x component of spins
extern std::vector<double> y; //saves the y component of spins
extern std::vector<double> z; //saves the z component of spins


extern std::vector<double> nx; //saves the x component of the uni-axial vector for all the spins
extern std::vector<double> ny; //saves the y component of the uni-axial vector for all the spins
extern std::vector<double> nz; //saves the z component of the uni-axial vector for all the spins



extern std::vector<double> k;  //saves the anisotropy constant for each spin


extern std::vector<double> lambda; //saves the lagrangian constraints

extern std::vector<double> global; //saves the global magnetisation


//namespace functions


void print_f();//function to print the spins' components

}





//namespace for structure related variables and functions
namespace structure {

//namespace variables



//namespace functions

void initial_structure_f(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z);//function to generate initial structure

void uniaxial_anisotropy_f(std::vector<double>& nx, std::vector<double>& ny, std::vector<double>& nz,
                           std::vector<double>& k); //function to create the uniaxial anisotropy of the system


void init_lagrangian_constraints_f(std::vector<double>& lambda); //function to initialise the lagrangian constraints

}

#endif // STRUCTURE_H
