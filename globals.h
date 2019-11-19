//GLOBALS MODULE - HEADER FILE//
///////////////////////////////


//----------------------------------- THIS HEADER FILE CONTAINS GLOBAL CONSTANTS -----------------------------------------//

#ifndef GLOBALS_H
#define GLOBALS_H


//C++ libraries
#include<cmath>



const int n1 = 100 ; //number of lines = number of atoms on a column   ###### for the 2d
const int n2 = 1 ; //number of columns  =  number of atoms on a line ######     case
const int n3 = 1 ; //number of atoms in depth

const int no_of_atoms = n1*n2*n3;


const int no_of_dimensions = 1;


const double J = 1.0;


const int uniform_init_structure = 1; // 1 = uniform structure, 0 = half & half
const int uniform_init_anisotropy = 1;// 1 = uniform anisotropy, 0 = half & half
const int uniform_lagrangian_constraints = 1;// 1 = uniform lagr. constraints -1 = automatically set the constraints
const int uniform_exchange_interaction = 1;// 1 = uniform exchange 0 = half&half

//----------ORIENTATION CONSTRAINT-------------//
const double v0_constrain_x=0.0;
const double v0_constrain_y=0.0;
const double v0_constrain_z=0.0;
//---------------------------------------------//


//--------------LAGRANGIAN PARAMETERS----------//
const double lambda_x = 1e-5;
const double lambda_y = 1e-5;
const double lambda_z = 1e-5;

const double l1 = 1e-5;

//---------------------------------------------//


//--------------FIELD COMPONENTS---------------//
const double hx = 0.0;
const double hy = 0.0;
const double hz = 1.0;

const double H = 0;

//---------------------------------------------//


//------------CONSTANTS------------------------//
const double mu_s = 4.4*1e-3; //spin magnetic moment


// Fe: 2.92*1e-3
// FePt: 4.4*1e-3
//---------------------------------------------//

#endif // GLOBALS_H
