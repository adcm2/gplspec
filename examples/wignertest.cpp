#include <iomanip>
#include <iostream>
#include <cmath>
#include <numbers>
#include <GSHTrans/All>
#include <Eigen/Dense>

int main (){
    using namespace GSHTrans;
    int lmax = 2;auto beta = 0.0;
    std::cout << "Type in lmax: ";
   std::cin >> lmax;
    std::cout << "Type in beta: ";
   std::cin >> beta;
    
    std::vector<Eigen::MatrixXcd> vec_wig;

    auto wigtemp = GSHTrans::Wigner<double, Ortho, All, All, Single, ColumnMajor>(lmax, lmax, lmax, beta);
    // std::vector<std::vector<std::vector<std::complex<double>>>> vec_wig;
    for (int l = 0; l < lmax + 1; ++l) {
        std::cout << "Working on l: " << l << "\n";
        auto normfactor = std::sqrt((4.0 * 3.1415926535897932) / (2 * l + 1));

        Eigen::MatrixXcd mat_tmp = Eigen::MatrixXcd::Zero(2 * l + 1, 2 * l + 1);
    //   std::vector<std::vector<std::complex<double>>> vec_tmp;
      // fill out matrix
      for (int N = -l; N < l + 1; ++N) {
        // std::vector<std::complex<double>> vec_tmp2;
         auto dl = wigtemp[N];
        //  auto dl = wigtemp[N];
         for (int m = -l; m < l + 1; ++m) {
            // std::cout << "l,N,m: " << l << "," << N << "," << m << "\n";
            auto tmpval = dl[l,m] * normfactor;
            mat_tmp(m + l, N+ l) = tmpval;
            // std::cout << std::setprecision(16) << tmpval<< "\n";
            //
                                    
                                    // vec_tmp2.push_back(tmpval);
           
         }
            // vec_tmp.push_back(vec_tmp2);
       
      }
        vec_wig.push_back(mat_tmp);
        std::cout << "Matrix for l=" << l << ":\n" << mat_tmp << "\n";
      std::cout << "\n\n";
        // vec_wig.push_back(vec_tmp);
   
   }
    return 0;}