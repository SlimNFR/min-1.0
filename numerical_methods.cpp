//NUM.METHODS MODULE - SOURCE FILE//
////////////////////////////////////


//C++ libraries
#include <vector>
#include <math.h>
#include <iostream>

//Code defined libraries
#include "globals.h"
#include "energy.h"
#include "utils.h"
#include "numerical_methods.h"
#include "structure.h"
#include "gradient.h"





void update_variable(std::vector<double>&x, std::vector<double>&y, std::vector<double>&z,std::vector<double>&lambda, std::vector<double> update_vec,double alpha)
{



for(int i =0; i<no_of_atoms;i++)
{

x[i] = x[i] - alpha * update_vec[i];
y[i] = y[i] - alpha * update_vec[i+no_of_atoms];
z[i] = z[i] - alpha * update_vec[i+2*no_of_atoms];

x[i] = x[i]/utils::compute_spin_magnitude(x[i],y[i],z[i]);
y[i] = y[i]/utils::compute_spin_magnitude(x[i],y[i],z[i]);
z[i] = z[i]/utils::compute_spin_magnitude(x[i],y[i],z[i]);


}


for(int i = 0; i<4 ; i++)
{


    lambda[i] = lambda[i] - alpha * update_vec[i+3*no_of_atoms];
}



}


//double steepest_descent (double function,std::vector<double>current_vars,std::vector<double>gradient, double tol, int max_iter)
//{

//    std::vector<double> temp_vars;
//    temp_vars.resize(4*no_of_atoms + 3);
//    double temp_function;


//    double alpha = 1e-3;//cat de mult mergi in directia data de gradient

//    k=0;

//    do{

//        for(unsigned int i = 0; i < current_vars.size();i++)
//        {

//            temp_vars [i] = current_vars[i] - alpha * gradient[i];

//        }



//    }while(function - temp_function <= tol && k >= max_iter);

//    while()
//    {
//        // aici trebuie sa chemi functia total_gradient



//        //aici faci current_vars = temp_vars
//        //calculeaza energia totala
//        //incrementeaza k


//    }






//}
