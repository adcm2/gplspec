#include <Eigen/Sparse>
#include <GaussQuad/All>
#include <Gravitational_Field/Test>
#include <Gravitational_Field/TestEllipticity>
#include <Gravitational_Field/Timer>
#include <chrono>
#include <fstream>
#include <iostream>
#include <string>

int
main() {
   using namespace GravityFunctions;
   // path to file:
   //    std::fstream prem;
   std::string pathtofile = "modeldata/prem.200";
   // std::string pathtofile = "modeldata/PREM500";

   Timer timer1;
   // timer1.start();
   auto parametermodel = EarthModels::EarthConstants<double>();
   TestTools::EarthModel myprem(pathtofile);
   // timer1.stop("Time to initialise model");

   // finding the gravitational acceleration:
   //    auto fourpig = 4.0 * 3.1415926535 * myprem.GravitationalConstant() *
   //                   myprem.DensityNorm();
   // timer1.start();

   //
   auto pi_db = 3.1415926535;
   auto bigGnorm = std::pow(myprem.LengthNorm(), 3.0) /
                   (myprem.MassNorm() * std::pow(myprem.TimeNorm(), 2.0));
   auto bigg = 6.6743 * std::pow(10.0, -11.0) / bigGnorm;
   auto eightpig = 8.0 * pi_db * bigg;
   auto fourpig = 4.0 * pi_db * bigg;

   // std::cout << myprem.GravitationalConstant() << "\n";
   // find the nodes
   double rad_bigstep = 200.0 / myprem.LengthNorm();
   double maxballrad = myprem.OuterRadius();
   std::vector<double> vec_noderadii =
       Radial_Node(myprem, rad_bigstep, maxballrad);
   std::vector<double> vec_hw(vec_noderadii.size() - 1);
   std::generate(vec_hw.begin(), vec_hw.end(),
                 [q = 0, &vec_noderadii]() mutable {
                    return vec_noderadii[q + 1] - vec_noderadii[q++];
                 });   // node widths

   ////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////

   // gravity
   std::vector<double> vec_gravity(vec_noderadii.size(), 0.0);
   auto itbeg = vec_gravity.begin();
   std::generate(
       vec_gravity.begin() + 1, vec_gravity.end(),
       [q = 0, laynum = 0, &myprem, &vec_hw, &vec_noderadii, itbeg]() mutable {
          if (vec_noderadii[q + 1] > myprem.UpperRadius(laynum)) {
             ++laynum;
          };
          return *(itbeg + q) +
                 0.5 * vec_hw[q] *
                     (myprem.Density(laynum)(vec_noderadii[q]) *
                          vec_noderadii[q] * vec_noderadii[q] +
                      myprem.Density(laynum)(vec_noderadii[++q]) *
                          vec_noderadii[q] * vec_noderadii[q]);
       });

   // generate Gauss grid for more accurate
   int npoly = 5;
   auto q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(npoly + 1);
   auto vec_allradii = All_Node(vec_noderadii, q);
   std::vector<double> vec_gravity_exact(vec_noderadii.size(), 0.0);
   // for (int idx = 0; idx < npoly + 1; ++idx) {
   //    std::cout << q.X(idx) << "\n";
   // }
   {
      auto laynum = 0;
      for (int idx = 1; idx < vec_noderadii.size(); ++idx) {
         if (vec_noderadii[idx] > myprem.UpperRadius(laynum)) {
            ++laynum;
         };
         vec_gravity_exact[idx] = vec_gravity_exact[idx - 1];
         for (int idxpoly = 0; idxpoly < npoly + 1; ++idxpoly) {
            auto newrad =
                0.5 *
                ((vec_noderadii[idx - 1] + vec_noderadii[idx]) +
                 q.X(idxpoly) * (vec_noderadii[idx] - vec_noderadii[idx - 1]));
            vec_gravity_exact[idx] +=
                q.W(idxpoly) * myprem.Density(laynum)(newrad) * newrad *
                newrad * (vec_noderadii[idx] - vec_noderadii[idx - 1]) / 2.0;
         }
      }
   }
   // get Earth mass
   auto earthM =
       4.0 * 3.1415926535 * vec_gravity_exact[vec_noderadii.size() - 1];
   std::cout << "Mass: " << earthM * myprem.MassNorm() << "\n";

   // divide by radius squared
   for (int idx = 1; idx < vec_noderadii.size(); ++idx) {
      vec_gravity[idx] *= fourpig / (vec_noderadii[idx] * vec_noderadii[idx]);
      vec_gravity_exact[idx] *=
          fourpig / (vec_noderadii[idx] * vec_noderadii[idx]);

      // vec_gravity[idx] *= 1.0 / (vec_noderadii[idx] * vec_noderadii[idx]);
      // vec_gravity_exact[idx] *= 1.0 / (vec_noderadii[idx] *
      // vec_noderadii[idx]);
   }
   std::cout << "Gravity at surface: " << vec_gravity_exact.back() << " "
             << vec_gravity_exact.back() * myprem.AccelerationNorm() << "\n";
   // timer1.stop("Time to find gravity");

   // ellipticity via forward difference solver
   using T = Eigen::Triplet<double>;
   using SpMat = Eigen::SparseMatrix<double>;
   std::vector<T> tripletlist;   // list of non-zeros coefficients the
                                 // constraints final point

   // debug
   // middle coefficients

   // first row
   tripletlist.push_back(T(0, 0, 1.0));
   tripletlist.push_back(T(0, 1, -1.0));
   {
      int laynum = 0;
      for (int idx = 1; idx < vec_noderadii.size() - 1; ++idx) {
         if (vec_noderadii[idx] > myprem.UpperRadius(laynum)) {
            ++laynum;
         };
         double valleft, valmid, valright;
         double anm, an0, anp;
         double bnm, bn0, bnp;

         // first order derivative coefficients
         auto multa = 1.0 / (vec_hw[idx - 1] + vec_hw[idx]);
         anm = -multa * vec_hw[idx] / vec_hw[idx - 1];
         an0 =
             (vec_hw[idx] - vec_hw[idx - 1]) / (vec_hw[idx - 1] * vec_hw[idx]);
         anp = multa * vec_hw[idx - 1] / vec_hw[idx];

         // second order derivative coefficients
         auto multb = 2.0 / (vec_hw[idx - 1] + vec_hw[idx]);
         bnm = multb / vec_hw[idx - 1];
         bnp = multb / vec_hw[idx];
         bn0 = -2.0 / (vec_hw[idx - 1] * vec_hw[idx]);

         // multiplication factor
         auto multmat = eightpig * myprem.Density(laynum)(vec_noderadii[idx]) /
                        vec_gravity_exact[idx];

         // values in matrix
         auto cnm = bnm + multmat * anm;
         auto cn0 = bn0;
         cn0 += multmat * (an0 + 1.0 / vec_noderadii[idx]);
         cn0 -= 6.0 / std::pow(vec_noderadii[idx], 2.0);
         auto cnp = bnp + multmat * anp;

         // adding to tripletlist
         tripletlist.push_back(T(idx, idx - 1, cnm));
         tripletlist.push_back(T(idx, idx, cn0));
         tripletlist.push_back(T(idx, idx + 1, cnp));
      }
   }

   // end points
   tripletlist.push_back(T(vec_noderadii.size() - 1, vec_noderadii.size() - 2,
                           -1.0 / vec_hw[vec_hw.size() - 1]));
   tripletlist.push_back(
       T(vec_noderadii.size() - 1, vec_noderadii.size() - 1,
         2.0 / (myprem.OuterRadius()) + 1.0 / vec_hw[vec_hw.size() - 1]));

   // set sparse matrix
   SpMat A(vec_noderadii.size(), vec_noderadii.size());
   A.setFromTriplets(tripletlist.begin(), tripletlist.end());
   A.makeCompressed();

   // rhs vector
   Eigen::VectorXd vec_rhs = Eigen::VectorXd::Zero(vec_noderadii.size());
   auto frequencynorm = 1.0 / myprem.TimeNorm();
   auto Omega = 2.0 * pi_db / (24.0 * 3600.0);
   Omega *= 1.0 / frequencynorm;
   auto earthrad = myprem.LengthNorm() * myprem.OuterRadius();

   // auto earthM = 5.972 * std::pow(10.0, 24.0);
   vec_rhs(vec_noderadii.size() - 1) =
       2.5 * std::pow(Omega * myprem.OuterRadius(), 2.0) / (bigg * earthM);

   // solve
   Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
   solver.compute(A);
   Eigen::VectorXd vec_ell = solver.solve(vec_rhs);
   ////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////
   // using spectral element method

   // setup new parameters
   double rad_step_spec = 30000.0 / myprem.LengthNorm();
   double maxballrad_spec = myprem.OuterRadius();
   std::vector<double> vec_noderadii_spec =
       Radial_Node(myprem, rad_step_spec, maxballrad_spec);
   auto vec_allradii_spec = All_Node(vec_noderadii_spec, q);

   std::vector<int> vec_num, vec_idx;
   for (int idx = 0; idx < myprem.NumberOfLayers(); ++idx) {
      double laydepth = myprem.UpperRadius(idx) - myprem.LowerRadius(idx);
      int numelements = std::ceil(laydepth / rad_step_spec);
      vec_num.push_back(numelements);
   }
   vec_idx.resize(myprem.NumberOfLayers());
   vec_idx[0] = vec_num[0] * npoly;
   for (int idx = 1; idx < myprem.NumberOfLayers(); ++idx) {
      vec_idx[idx] = vec_idx[idx - 1] + vec_num[idx] * npoly;
   }

   // check
   // {
   //    auto idx1 = 0;
   //    for (auto &idx : vec_idx) {
   //       std::cout << idx << " " << vec_allradii_spec[idx] << " "
   //                 << myprem.LowerRadius(idx1) << " "
   //                 << myprem.UpperRadius(idx1++) << "\n";
   //    }
   // }

   // finding derivative of i-th lagrange interpolant at j-th point
   Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> _mat_gaussderiv;
   _mat_gaussderiv.resize(npoly + 1, npoly + 1);
   {
      auto pleg = Interpolation::LagrangePolynomial(q.Points().begin(),
                                                    q.Points().end());
      for (int idxi = 0; idxi < npoly + 1; ++idxi) {
         for (int idxj = 0; idxj < npoly + 1; ++idxj) {
            _mat_gaussderiv(idxi, idxj) = pleg.Derivative(idxi, q.X(idxj));
         }
      }
   }

   // finding gravity
   std::vector<double> vec_gravity_spec(vec_allradii_spec.size(), 0.0);

   {
      auto laynum = 0;
      for (int idx = 1; idx < vec_allradii_spec.size(); ++idx) {

         vec_gravity_spec[idx] = vec_gravity_spec[idx - 1];
         for (int idxpoly = 0; idxpoly < npoly + 1; ++idxpoly) {
            auto newrad =
                0.5 * ((vec_allradii_spec[idx - 1] + vec_allradii_spec[idx]) +
                       q.X(idxpoly) * (vec_allradii_spec[idx] -
                                       vec_allradii_spec[idx - 1]));
            vec_gravity_spec[idx] +=
                q.W(idxpoly) * myprem.Density(laynum)(newrad) * newrad *
                newrad * (vec_allradii_spec[idx] - vec_allradii_spec[idx - 1]) /
                2.0;
         }
         if (idx == vec_idx[laynum]) {
            // std::cout << vec_allradii_spec[idx] << " "
            //           << myprem.UpperRadius(laynum) << "\n";
            ++laynum;
         };
      }
   }

   // get Earth mass
   auto earthM_spec =
       4.0 * 3.1415926535 * vec_gravity_spec[vec_allradii_spec.size() - 1];
   std::cout << "Mass: " << earthM_spec << "\n";

   // divide by radius squared
   for (int idx = 1; idx < vec_allradii_spec.size(); ++idx) {
      auto radval = vec_allradii_spec[idx] * vec_allradii_spec[idx];
      vec_gravity_spec[idx] *= fourpig / radval;
   }

   // check gravity
   std::cout << "Gravity: "
             << vec_gravity_spec.back() * myprem.AccelerationNorm() << "\n";

   int nelem_spec = vec_noderadii_spec.size() - 1;   // total number of elements
   std::vector<double> nodespacing(nelem_spec), invnodespacingtwo(nelem_spec);
   for (int idxelem = 0; idxelem < nelem_spec; ++idxelem) {
      nodespacing[idxelem] =
          vec_noderadii_spec[idxelem + 1] - vec_noderadii_spec[idxelem];
      invnodespacingtwo[idxelem] = 2.0 / nodespacing[idxelem];
   };

   // finding value of density radially and then derivative
   std::vector<double> vec_density_spec(nelem_spec * (npoly + 1), 0.0);
   {
      auto laynum = 0;
      for (int idxelem = 0; idxelem < nelem_spec; ++idxelem) {
         for (int idxpoly = 0; idxpoly < npoly + 1; ++idxpoly) {
            // global index
            auto idxglobal = idxelem * (npoly + 1) + idxpoly;
            auto crad = vec_allradii_spec[idxelem * npoly + idxpoly];
            vec_density_spec[idxglobal] = myprem.Density(laynum)(crad);

            if ((idxelem * npoly + idxpoly) == vec_idx[laynum]) {
               ++laynum;
            };
         }
      }
   }
   std::cout << "Density: "
             << vec_density_spec[nelem_spec * (npoly + 1) - 1] *
                    myprem.DensityNorm()
             << " "
             << vec_density_spec[nelem_spec * (npoly + 1) - 1] *
                    myprem.DensityNorm()
             << "\n";

   // finding radial derivative of density
   std::vector<double> vec_deriv_density(nelem_spec * (npoly + 1), 0.0);
   for (int idxelem = 0; idxelem < nelem_spec; ++idxelem) {
      for (int idxpoly = 0; idxpoly < npoly + 1; ++idxpoly) {
         // global index
         auto idxglobal = idxelem * (npoly + 1) + idxpoly;
         for (int idxk = 0; idxk < npoly + 1; ++idxk) {
            // sum of the coefficient of the k-th interpolating polynomial at
            // the idxpoly-th point
            auto idxgk = idxelem * (npoly + 1) + idxk;
            vec_deriv_density[idxglobal] +=
                vec_density_spec[idxgk] * _mat_gaussderiv(idxk, idxpoly);
            if (idxelem == nelem_spec - 1 && idxpoly == npoly) {
               // std::cout << vec_density_spec[idxgk] << "\n";
            }
         }
         vec_deriv_density[idxglobal] *= invnodespacingtwo[idxelem];
      }
   }
   std::cout << "Density derivative: "
             << vec_deriv_density.back() * myprem.DensityNorm() /
                    myprem.LengthNorm()
             << "\n";

   // using T = Eigen::Triplet<CFLOAT>;
   std::vector<T> spec_tripletlist;
   for (int idxelem = 0; idxelem < nelem_spec; ++idxelem) {
      // need a global index pointing to the "0-th" element in the box
      int idxglobal = idxelem * npoly;
      for (int idxi = 0; idxi < npoly + 1; ++idxi) {
         for (int idxj = 0; idxj < npoly + 1; ++idxj) {
            // finding (i,j) element
            auto tmp = 0.0;
            for (int idxk = 0; idxk < npoly + 1; ++idxk) {
               tmp += q.W(idxk) *
                      std::pow(vec_allradii_spec[idxglobal + idxk], 2.0) *
                      _mat_gaussderiv(idxi, idxk) *
                      _mat_gaussderiv(idxj, idxk) * invnodespacingtwo[idxelem];
            }
            // tmp *= invnodespacingtwo[idxelem];

            // adding to list
            spec_tripletlist.push_back(
                T(idxglobal + idxi, idxglobal + idxj, tmp));
         }

         // diagonal component
         int idxdiag = idxglobal + idxi;
         spec_tripletlist.push_back(
             T(idxdiag, idxdiag, 6.0 * 0.5 * nodespacing[idxelem] * q.W(idxi)));

         // term 2
         if (idxdiag != 0) {
            int idx_layers = idxelem * (npoly + 1) + idxi;
            auto tmp = 0.0;
            tmp += std::pow(vec_allradii_spec[idxdiag], 2.0);
            tmp *= vec_deriv_density[idx_layers];
            tmp *= 1.0 / vec_gravity_spec[idxdiag];
            tmp *= fourpig * nodespacing[idxelem] * 0.5;
            tmp *= q.W(idxi);
            spec_tripletlist.push_back(T(idxdiag, idxdiag, tmp));
         }
      }
   }

   // add in external boundary :
   {
      int idxt = nelem_spec * npoly;
      spec_tripletlist.push_back(
          T(idxt, idxt,
            myprem.OuterRadius() *
                (3.0 - fourpig * myprem.OuterRadius() *
                           vec_density_spec.back() / vec_gravity_spec.back())));
   }

   // add in internal boundaries
   {
      for (int idx = 0; idx < myprem.NumberOfLayers() - 1; ++idx) {

         auto tmp = myprem.Density(idx + 1)(myprem.LowerRadius(idx + 1)) -
                    myprem.Density(idx)(myprem.UpperRadius(idx));
         tmp *= fourpig;
         // tmp *= 2.0;
         tmp *= 1.0 / vec_gravity_spec[vec_idx[idx]];
         tmp *= std::pow(vec_allradii_spec[vec_idx[idx]], 2.0);
         spec_tripletlist.push_back(T(vec_idx[idx], vec_idx[idx], tmp));
      }
   }

   Eigen::SparseMatrix<double> mat_mass(vec_allradii_spec.size(),
                                        vec_allradii_spec.size());
   // std::vector<T> vec_settest;
   // vec_settest.push_back(T(0, 0, 0.1));
   // std::cout << "Hello\n";
   // mat_mass.setFromTriplets(vec_settest.begin(), vec_settest.end());
   mat_mass.setFromTriplets(spec_tripletlist.begin(), spec_tripletlist.end());
   // std::cout << "Hello\n";
   ////////////////////////////////////////////////////////////////////
   // finding force

   Eigen::VectorXd vec_force = Eigen::VectorXd::Zero(vec_allradii_spec.size());

   auto prefactor = fourpig * Omega * Omega * sqrt(4.0 * pi_db / 5.0) / 3.0;
   {

      int laynum = 0;
      // int idx_glob = 1;
      for (int idxelem = 0; idxelem < nelem_spec; ++idxelem) {
         int idxpolymin = 0;
         // int idxpolymax = npoly;

         if (idxelem == 0) {
            idxpolymin = 1;
         }

         for (int idxpoly = idxpolymin; idxpoly < npoly + 1; ++idxpoly) {
            int idx_glob = idxelem * npoly + idxpoly;
            int idx_layers = idxelem * (npoly + 1) + idxpoly;
            auto crad = vec_allradii_spec[idx_glob];
            auto tmp = std::pow(vec_allradii_spec[idx_glob], 4.0);
            tmp *= vec_deriv_density[idx_layers] / vec_gravity_spec[idx_glob];
            tmp *= -prefactor * nodespacing[idxelem] / 2.0 * q.W(idxpoly);
            vec_force(idx_glob) += tmp;
            if ((idxelem * npoly + idxpoly) == vec_idx[laynum]) {
               ++laynum;
            };
         }
      }
   }

   // add in external boundary:
   {
      int idxt = nelem_spec * npoly;
      auto tmp = std::pow(myprem.OuterRadius(), 4.0);
      tmp *= vec_density_spec.back() / vec_gravity.back();
      tmp *= prefactor;
      vec_force(idxt) += tmp;
   }

   // add in internal boundaries
   {
      for (int idx = 0; idx < myprem.NumberOfLayers() - 1; ++idx) {
         // auto drho = myprem.Density(idx + 1)(myprem.LowerRadius(idx + 1)) -
         //             myprem.Density(idx)(myprem.UpperRadius(idx));
         auto tmp = myprem.Density(idx + 1)(myprem.LowerRadius(idx + 1)) -
                    myprem.Density(idx)(myprem.UpperRadius(idx));
         tmp *= -prefactor;
         tmp *= 1.0 / vec_gravity_spec[vec_idx[idx]];
         tmp *= std::pow(vec_allradii_spec[vec_idx[idx]], 4.0);
         vec_force(vec_idx[idx]) += tmp;
      }
   }

   Eigen::VectorXd vec_force2 = Eigen::VectorXd::Zero(vec_allradii_spec.size());
   vec_force2(vec_allradii_spec.size() - 1) =
       5.0 / 3.0 * sqrt(4.0 * pi_db / 5.0) * Omega * Omega *
       std::pow(myprem.OuterRadius(), 3.0);

   Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>, Eigen::Lower,
                         Eigen::AMDOrdering<int>>
       chol_solver;
   Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
       lu_solver;
   mat_mass.makeCompressed();
   chol_solver.compute(mat_mass);
   lu_solver.compute(mat_mass);
   Eigen::VectorXd vec_sol = lu_solver.solve(vec_force);
   Eigen::VectorXd vec_sol2 = lu_solver.solve(vec_force2);

   // forgot particular solution!!
   {
      auto prefact2 = Omega * Omega * sqrt(4.0 * pi_db / 5.0) / 3.0;
      for (int idx = 1; idx < vec_allradii_spec.size(); ++idx) {
         vec_sol(idx) +=
             prefact2 * vec_allradii_spec[idx] * vec_allradii_spec[idx];
      }
   }

   // ellipticity
   Eigen::VectorXd vec_ell2 = Eigen::VectorXd::Zero(vec_allradii_spec.size());
   for (int idx = 1; idx < vec_allradii_spec.size(); ++idx) {
      auto tmp = sqrt(4.0 * pi_db / 5.0) * 2.0 / 3.0;
      tmp *= vec_allradii_spec[idx] * vec_gravity_spec[idx];
      vec_ell2(idx) = vec_sol2(idx) / tmp;
   }
   // ellipticity at zero:
   for (int idxpoly = 1; idxpoly < npoly + 1; ++idxpoly) {
      vec_ell2(0) -= _mat_gaussderiv(idxpoly, 0) / _mat_gaussderiv(0, 0) *
                     vec_ell2(idxpoly);
   }

   ////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////
   // spectral element on Clairaut's equation
   std::vector<T> clairaut_tripletlist;
   for (int idxelem = 0; idxelem < nelem_spec; ++idxelem) {
      // need a global index pointing to the "0-th" element in the box
      int idxglobal = idxelem * npoly;
      for (int idxi = 0; idxi < npoly + 1; ++idxi) {
         int idxgi = idxglobal + idxi;
         int idx_density = idxelem * (npoly + 1) + idxi;
         for (int idxj = 0; idxj < npoly + 1; ++idxj) {
            // finding (i,j) element
            auto tmp = 0.0;
            for (int idxk = 0; idxk < npoly + 1; ++idxk) {
               tmp -= q.W(idxk) *
                      std::pow(vec_allradii_spec[idxglobal + idxk], 2.0) *
                      _mat_gaussderiv(idxi, idxk) * _mat_gaussderiv(idxj, idxk);
            }
            tmp *= invnodespacingtwo[idxelem];

            // adding to list
            clairaut_tripletlist.push_back(T(idxgi, idxglobal + idxj, tmp));

            // add in kappa * \dot{\epsilon} term
            if (idxgi != 0) {
               clairaut_tripletlist.push_back(T(
                   idxgi, idxglobal + idxj,
                   2.0 * fourpig * std::pow(vec_allradii_spec[idxgi], 2.0) *
                       vec_density_spec[idx_density] / vec_gravity_spec[idxgi] *
                       _mat_gaussderiv(idxj, idxi) * q.W(idxi)));
               clairaut_tripletlist.push_back(
                   T(idxgi, idxglobal + idxj,
                     -2.0 * vec_allradii_spec[idxgi] *
                         _mat_gaussderiv(idxj, idxi) * q.W(idxi)));
            }
         }

         // diagonal component
         // int idxdiag = idxglobal + idxi;
         if (idxgi != 0) {
            clairaut_tripletlist.push_back(
                T(idxgi, idxgi,
                  0.5 * nodespacing[idxelem] * 2.0 * fourpig *
                      vec_allradii_spec[idxgi] * vec_density_spec[idx_density] /
                      vec_gravity_spec[idxgi] * q.W(idxi)));
            clairaut_tripletlist.push_back(
                T(idxgi, idxgi, -6.0 * q.W(idxi) * 0.5 * nodespacing[idxelem]));
         }
      }
   }

   // add in final point
   clairaut_tripletlist.push_back(T(vec_allradii_spec.size() - 1,
                                    vec_allradii_spec.size() - 1,
                                    -2.0 * myprem.OuterRadius()));

   //////////////////////////////////////////////////////
   // second way for Clairaut
   //  test
   std::vector<T> clairaut_tripletlist2;
   for (int idxelem = 0; idxelem < nelem_spec; ++idxelem) {
      // need a global index pointing to the "0-th" element in the box
      int idxglobal = idxelem * npoly;
      for (int idxi = 0; idxi < npoly + 1; ++idxi) {
         int idxgi = idxglobal + idxi;
         int idx_density = idxelem * (npoly + 1) + idxi;
         for (int idxj = 0; idxj < npoly + 1; ++idxj) {
            // finding (i,j) element
            // if (idxgi != 0) {
            auto tmp = 0.0;
            for (int idxk = 0; idxk < npoly + 1; ++idxk) {
               tmp -= q.W(idxk) * _mat_gaussderiv(idxi, idxk) *
                      _mat_gaussderiv(idxj, idxk);
            }
            tmp *= invnodespacingtwo[idxelem];

            // adding to list
            clairaut_tripletlist2.push_back(T(idxgi, idxglobal + idxj, tmp));
            // }

            // }
            if (idxgi == 0) {
               for (int idxk = 0; idxk < npoly + 1; ++idxk) {
                  // clairaut_tripletlist2.push_back(
                  //     T(0, idxk, _mat_gaussderiv(idxk, 0)));
               }

               // clairaut_tripletlist2.push_back(T(0, 1, -1.0));
            }
            // add in kappa * \dot{\epsilon} term
            if (idxgi != 0) {
               clairaut_tripletlist2.push_back(
                   T(idxgi, idxglobal + idxj,
                     2.0 * fourpig * vec_density_spec[idx_density] /
                         vec_gravity_spec[idxgi] * _mat_gaussderiv(idxj, idxi) *
                         q.W(idxi)));
               // clairaut_tripletlist2.push_back(
               //     T(idxgi, idxglobal + idxj,
               //       -2.0 / vec_allradii_spec[idxgi] *
               //           _mat_gaussderiv(idxj, idxi) * q.W(idxi)));
            }
         }

         // diagonal component
         if (idxgi != 0) {
            clairaut_tripletlist2.push_back(
                T(idxgi, idxgi,
                  0.5 * nodespacing[idxelem] * 2.0 * fourpig /
                      vec_allradii_spec[idxgi] * vec_density_spec[idx_density] /
                      vec_gravity_spec[idxgi] * q.W(idxi)));
            clairaut_tripletlist2.push_back(
                T(idxgi, idxgi,
                  -6.0 * q.W(idxi) * 0.5 * nodespacing[idxelem] /
                      std::pow(vec_allradii_spec[idxgi], 2.0)));
         }
      }
   }

   // add in final point
   clairaut_tripletlist2.push_back(T(vec_allradii_spec.size() - 1,
                                     vec_allradii_spec.size() - 1,
                                     -2.0 / myprem.OuterRadius()));

   Eigen::SparseMatrix<double> mat_cl(vec_allradii_spec.size(),
                                      vec_allradii_spec.size()),
       mat_cl2(vec_allradii_spec.size(), vec_allradii_spec.size());
   mat_cl.setFromTriplets(clairaut_tripletlist.begin(),
                          clairaut_tripletlist.end());
   mat_cl.makeCompressed();
   mat_cl2.setFromTriplets(clairaut_tripletlist2.begin(),
                           clairaut_tripletlist2.end());
   mat_cl2.makeCompressed();
   std::cout << "0,0: " << mat_cl.coeff(0, 0) << "\n";

   Eigen::VectorXd vec_fcl = Eigen::VectorXd::Zero(vec_allradii_spec.size());
   vec_fcl(vec_allradii_spec.size() - 1) = -5.0 / 2.0 * Omega * Omega *
                                           std::pow(myprem.OuterRadius(), 4.0) /
                                           (bigg * earthM_spec);
   Eigen::VectorXd vec_fcl2 = Eigen::VectorXd::Zero(vec_allradii_spec.size());
   vec_fcl2(vec_allradii_spec.size() - 1) =
       -5.0 / 2.0 * Omega * Omega * std::pow(myprem.OuterRadius(), 2.0) /
       (bigg * earthM_spec);

   Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
       lusolver, lusolver2;

   lusolver.compute(mat_cl);
   Eigen::VectorXd vec_sol_cl = lusolver.solve(vec_fcl);

   lusolver2.compute(mat_cl2);
   Eigen::VectorXd vec_sol_cl2 = lusolver2.solve(vec_fcl2);

   ////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////
   auto pathtooutputfile = "./work/gravitycheck.out";
   auto file2 = std::ofstream(pathtooutputfile, std::ios::out);

   for (int idx = 0; idx < vec_noderadii.size(); ++idx) {

      file2 << std::setprecision(16) << vec_noderadii[idx] * myprem.LengthNorm()
            << " " << vec_gravity[idx] * myprem.AccelerationNorm() << " "
            << vec_gravity_exact[idx] * myprem.AccelerationNorm() << " "
            << vec_ell(idx) << "\n";
   }
   file2.close();

   pathtooutputfile = "./work/ellipticity2.out";
   file2 = std::ofstream(pathtooutputfile, std::ios::out);
   int laynum = 0;
   for (int idxelem = 0; idxelem < nelem_spec; ++idxelem) {
      int idxpolymax = npoly + 1;
      if (idxelem == nelem_spec - 1) {
         idxpolymax = npoly + 1;
      }
      for (int idxpoly = 0; idxpoly < idxpolymax; ++idxpoly) {
         int idx = idxelem * npoly + idxpoly;
         file2 << std::setprecision(16)
               << vec_allradii_spec[idx] * myprem.LengthNorm() << " "
               << vec_gravity_spec[idx] * myprem.AccelerationNorm() << " "
               << vec_density_spec[idxelem * (npoly + 1) + idxpoly] *
                      myprem.DensityNorm()
               << " "
               << myprem.Density(laynum)(vec_allradii_spec[idx]) *
                      myprem.DensityNorm()
               << " " << vec_sol_cl2(idx) << " " << vec_sol_cl(idx) << " "
               << vec_ell2(idx) << "\n";
         if (idx == vec_idx[laynum]) {
            ++laynum;
         };
      }
   }

   return 0;
   // return 0;
};
