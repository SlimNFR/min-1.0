#ifndef NUMERICAL_METHODS_H
#define NUMERICAL_METHODS_H


void test_dir(double alpha);
void update_variable(std::vector<double>&x, std::vector<double>&y, std::vector<double>&z,
                     std::vector<double>&lambda,
                     std::vector<double> update_vec,double alpha); //function to update the variables




#endif // NUMERICAL_METHODS_H
