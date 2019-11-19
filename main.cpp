//MAIN MODULE - SOURCE FILE///
//////////////////////////////


//----------------------------------------------------------------APPEALING NAME FOR ME CODE: BLABLA V.1.0-------------------------------------------------------//


//C++ libraries
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

//Code defined libraries
#include "energy.h"
#include "globals.h"
#include "utils.h"
#include "structure.h"
#include "gradient.h"
#include "numerical_methods.h"

int main()
{

    std::ofstream f1;
    std::ofstream f2;


    f1.open("spins.txt", std::ofstream::out);
    f2.open("energy_and_gradient.txt",std::ofstream::out);


      //--------------------------------------- initialise problem --------------------//
      structure::initial_structure_f(spin::x,spin::y,spin::z); //generate initial spin structure
      //spin::print_f(); // print the spins' components

      structure::uniaxial_anisotropy_f(spin::nx,spin::ny,spin::nz,spin::k);  //generate the anisotropy configuration
      structure::init_lagrangian_constraints_f(spin::lambda); //setup the initial lagrangian constraints



      std::cout<<"SPINI: "<<spin::x[0]<<","<<spin::y[0]<<","<<spin::z[0]<<"\n";

//      utils::save_state(utils::temp_x,utils::temp_y,utils::temp_z,utils::temp_lambda,utils::temp_gradient,utils::temp_energy);



      //---------------------------------------compute total energy and total gradient --------------------//

//      energy::total_f();
//      gradient::total_f();
//      std::cout<<utils::compute_modul(gradient::total);




     // ---------------------------------------------------test steepest descent--------------------------------------//





            std::vector<double>s_vec;//the search vector
            s_vec.resize(3*no_of_atoms+4);
            double temp_en;


            double max_iter = 1e4;
            //double diff;
            double tol = -1;
            int k = 1;




            double alpha = 1e-1;




            do
            {




                energy::total_f();
                temp_en =energy::total;
                gradient::total_f();
                utils::compute_global_magnetisation();
                utils::print_all(f1,f2,k);



//std::cout<<spin::lambda[0]<<" "<<spin::lambda[1]<<" "<<spin::lambda[2]
                            for(unsigned int i = 0; i<3*no_of_atoms+4;i++)
                                {
                                 std::cout<<"grad.component"<<gradient::total[i]<<"\n";


                                s_vec[i] = gradient::total[i]/(utils::compute_modul(gradient::total));
}


//                                std::cout<<s_vec[i]<<" ";





                            update_variable(spin::x,spin::y,spin::z,spin::lambda,s_vec,alpha);
                            energy::total_f();


////                            while(energy::total > temp_en)
////                            {
////                                update_variable(spin::x,spin::y,spin::z,spin::lambda,s_vec,alpha);
////                                energy::total_f();


////                            }



//                               gradient::total_f();


////                               diff = energy::total - temp_en ;
////                               if(diff > 0)
////                                   alpha = alpha*1e1;



//                              /* //diff = energy::total_energy - temp_en;*/

                               k++;
                               std::cout<<"PASUL "<<k<<std::endl;
//                               std::cout<<diff<<","<<temp_en<<","<<energy::total<<"\n";
//                               std::cout<<"GRADIENT:r "<<gradient::lagrangian_x[0]<<" "<<gradient::lagrangian_y[0]<<" "<<gradient::lagrangian_z[0];





            }while(k <= max_iter && abs(energy::total - temp_en)>tol && abs(spin::global[2] - v0_constrain_z)>tol);




            f1.close();
            f2.close();

    return 0;
}
