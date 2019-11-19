//UTILS MODULE - HEADER FILE//
//////////////////////////////

#ifndef UTILS_H
#define UTILS_H


//C++ libraries
#include <vector>
#include <fstream>


//Code defined libraries
namespace utils {


//namespace variables

extern std::vector<double>global_magnetisation; // array to save the global magnetisation vector's components
extern std::vector<double> sums; //array to save the sums of the spins' components on x,y,z and the sum of the x components squared + sum of the y components squared + sum of the z components squared

extern std::vector<double>temp_x;
extern std::vector<double>temp_y;
extern std::vector<double>temp_z;
extern std::vector<double>temp_lambda;       //temp variables to store the current state
extern std::vector<double>temp_gradient;
extern double temp_energy;


//namespace functions


void print_all(std::ofstream& f1, std::ofstream& f2,int step);

void save_state(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z,
                std::vector<double>& lambda,
                std::vector<double>& gradient,
                double& energy);


double compute_modul(std::vector<double> vec);


void compute_global_magnetisation(); //function to compute the global magnetisation
void compute_magnetisation_sums(std::vector<double>x, std::vector<double> y, std::vector<double>z,
                                std::vector<double>& sums); //function to compute the sums of the spin components on x,y,z and the  (m_ix)^2 + (m_iy)^2 + (m_iz)^2  term

double compute_spin_magnitude(double x, double y, double z);



}


void compute_gradient(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z,
                      std::vector<double>& lambda,
                      std::vector<double>& gradient);




void create_init_lagrangian_constraints(std::vector<double>& lambda); //setup the initial lagrangian constraints


double compute_lagrangian_term(std::vector<double> x, std::vector<double> y, std::vector<double> z,
                               std::vector<double> lambda); //function to compute the lagrangian term
                                                                                                                                     //NOTE!!! the first three elements in the
                                                                                                                                     //lambda array, belong to the vector LAMBDA









void create_uniax_anis(std::vector<double>& nx, std::vector<double>& ny, std::vector<double>& nz,
                       std::vector<double>& k); //assign values to the nx,ny,nz,k vectors


void generate_init_structure(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z); //generate the initial spin config.


double compute_exchange_energy(std::vector<double> x, std::vector<double> y, std::vector<double> z); //function to calculate the exchange energy














#endif // UTILS_H
