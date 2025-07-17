#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <GaussQuad/All>
#include <Interpolation/All>
#include <PlanetaryModel/All>
// #include <TomographyModels/All>
// #include <FFTWpp/All>
// #include <GSHTrans/All>

#include <algorithm>
#include <complex>
#include <concepts>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <chrono>
#include <random>
#include <vector>
#include "Spherical_Integrator.h"

int main()
{
    using namespace std::chrono;

    using namespace Interpolation;
    using namespace EarthModels;

    // polynomial order
    int npoly = 6;
    std::cout << "Order of polynomial: \n";
    std::cin >> npoly;

    ///////////////////////////////////////////////////////////////

    auto myprem = PREM();

    // number of layers and scale factor
    int nlayer = myprem.NumberOfLayers();
    double rad_scale = myprem.OuterRadius();

    // number of nodes within each layer, including those at boundaries
    std::vector<int> vec_numnodes(nlayer, 6);
    vec_numnodes[1] = 6;

    // node values:
    auto start = high_resolution_clock::now();
    std::vector<double> vec_noderadii;

    // add first node at r = 0
    vec_noderadii.push_back(0.0);
    for (int idx = 0; idx < nlayer; ++idx)
    {
        double dx = (myprem.UpperRadius(idx) - myprem.LowerRadius(idx)) / (static_cast<double>(vec_numnodes[idx] - 1) * rad_scale);
        for (int idx2 = 1; idx2 < vec_numnodes[idx]; ++idx2)
        {

            vec_noderadii.push_back(myprem.LowerRadius(idx) / rad_scale + dx * static_cast<double>(idx2));
        }
    }

    int nelem = vec_noderadii.size() - 1; // total number of elements

    // int npoly = 6;                  // number of polynomials within each element
    int matlen = nelem * npoly + 1; // size of matrix

    // function:
    int N = 5;
    std::vector<double> x(N), y(N);

    // // finite element matrices
    Eigen::VectorXd vecforce = Eigen::VectorXd::Zero(matlen);
    Eigen::VectorXd vecsol = Eigen::VectorXd::Zero(matlen);
    Eigen::SparseMatrix<double> matspec(matlen, matlen);
    using T = Eigen::Triplet<double>;
    std::vector<T> tripletList;
    // std::vector<T> tripletList(nelem * (npoly + 1) * (npoly + 1) + 1);
    tripletList.reserve(nelem * (npoly + 1) * (npoly + 1) + 1);

    // generate Gauss grid
    auto q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(npoly + 1);

    // vector of Lagrange polynomials:
    auto pleg = LagrangePolynomial(q.Points().begin(), q.Points().end());

    auto posforce = [](double x)
    { return std::cos(x); };
    auto func = [&pleg](int i, int j, double xval)
    {
        return pleg(i, xval) * pleg(j, xval);
    };
    auto func1 = [&pleg](int i, int j, double xval)
    {
        return pleg.Derivative(i, xval) * pleg(j, xval);
    };
    auto func2 = [&pleg](int i, int j, double xval)
    {
        return pleg.Derivative(i, xval) * pleg.Derivative(j, xval);
    };
    auto rphys = [&vec_noderadii](int idxelem, double x)
    {
        return ((vec_noderadii[idxelem + 1] - vec_noderadii[idxelem]) * x + (vec_noderadii[idxelem + 1] + vec_noderadii[idxelem])) *
               0.5;
    };
    auto rscale = [&vec_noderadii](int idxelem, double x)
    {
        return (x + (vec_noderadii[idxelem + 1] + vec_noderadii[idxelem]) / (vec_noderadii[idxelem + 1] - vec_noderadii[idxelem]));
    };

    // l
    int lval = 0;
    const double zetasq = static_cast<double>(lval) * static_cast<double>(lval + 1);

    // initialising vector of matrices which will have the values of the derivatives of the Lagrange interpolant at the Gaussian quadrature points
    std::vector<Eigen::MatrixXd> vec_mat_lderiv;
    for (int idxk = 0; idxk < npoly + 1; ++idxk)
    {
        Eigen::MatrixXd mat_lderiv = Eigen::MatrixXd::Zero(npoly + 1, npoly + 1);
        for (int idxi = 0; idxi < npoly + 1; ++idxi)
        {
            for (int idxj = 0; idxj < npoly + 1; ++idxj)
            {
                mat_lderiv(idxi, idxj) = func2(idxi, idxj, q.X(idxk));
            };
        };
        vec_mat_lderiv.push_back(mat_lderiv);
    }

    // lambda to find value of r^2 l_i'(x_k) l_j'(x_k), with x_k as the Gaussian quadrature point
    auto funcnew = [&rphys, &rscale, &vec_mat_lderiv, &q](int idxelem, int idxi, int idxj, int idxk)
    {
        return rphys(idxelem, q.X(idxk)) * rscale(idxelem, q.X(idxk)) * vec_mat_lderiv[idxk](idxi, idxj);
    };
    double mytest;
    std::vector<double> intvec(npoly + 1);

    // finding the integrals using the inner product and the vector with the derivative values at the interpolation points
    for (int idxelem = 0; idxelem < nelem; ++idxelem)
    {
        for (int idxi = 0; idxi < npoly + 1; ++idxi)
        {
            for (int idxj = 0; idxj < npoly + 1; ++idxj)
            {
                // value of function at Gaussian quadrature points
                for (int idxpoly = 0; idxpoly < npoly + 1; ++idxpoly)
                {
                    intvec[idxpoly] = funcnew(idxelem, idxi, idxj, idxpoly);
                };

                // perform integration, different result if i = j or not
                if (idxi == idxj)
                {
                    tripletList.push_back(T(idxelem * npoly + idxi, idxelem * npoly + idxj, -std::inner_product(intvec.begin(), intvec.end(), q.Weights().cbegin(), 0.0) + zetasq * q.W(idxi)));
                }
                else
                {

                    tripletList.push_back(T(idxelem * npoly + idxi, idxelem * npoly + idxj, -std::inner_product(intvec.begin(), intvec.end(), q.Weights().cbegin(), 0.0)));
                }
            }
        }
    }

    // final point
    tripletList.push_back(T(matlen - 1, matlen - 1, -static_cast<double>(lval + 1) * vec_noderadii[nelem]));

    // constructing the sparse matrix and compressing it
    matspec.setFromTriplets(tripletList.begin(), tripletList.end());
    matspec.makeCompressed();
    // std::cout << "10, 8 member: " << matspec(10, 8) << std::endl;

    // std::cout << "Size of matrix: " << matlen << std::endl;

    // auto stop = high_resolution_clock::now();
    // auto duration = duration_cast<microseconds>(stop - start);
    // std::cout << "Setup time pt 1: "
    //           << duration.count() / 1000000.0 << " seconds" << std::endl;

    // Now to find the force vector
    // start = high_resolution_clock::now();

    // values of G, pi and the scales relevant
    const double bigg_db = 6.6743 * std::pow(10.0, -11.0);
    const double pi_db = 3.1415926535;
    double scdiff = rad_scale / myprem.OuterRadius();
    double multfact = 4.0 * pi_db * bigg_db * std::pow(rad_scale, 2.0); // multiplication factor
    int laynum = 0;                                                     // initialise layer number

    auto funcforforce = [&myprem, &rphys, &rscale, &scdiff, &q, &vec_noderadii](int idxelem, int idxpoly, int laynum)
    {
        return std::pow(0.5 * (vec_noderadii[idxelem + 1] - vec_noderadii[idxelem]), 3.0) * std::pow(rscale(idxelem, q.X(idxpoly)), 2.0) * myprem.Density(laynum)(rphys(idxelem, q.X(idxpoly)) * scdiff);
        // return rphys(idxelem, x) * rphys(idxelem, x) * pleg(idxpoly, x) * 1000.0; // test case with constant density
    };
    // force vector
    for (int idxelem = 0; idxelem < nelem; ++idxelem)
    {
        // finding the layer number
        if (!(vec_noderadii[idxelem] < myprem.UpperRadius(laynum) / rad_scale))
        {
            laynum += 1;
        };

        for (int idxpoly = 0; idxpoly < npoly + 1; ++idxpoly)
        {
            vecforce(idxelem * npoly + idxpoly) +=
                multfact * funcforforce(idxelem, idxpoly, laynum) * q.W(idxpoly);
        };
    };

    // stop = high_resolution_clock::now();
    // duration = duration_cast<microseconds>(stop - start);
    // std::cout << "Setup time pt 2: "
    //           << duration.count() / 1000000.0 << " seconds" << std::endl;

    ////////////////////////////////////////////////////////////////
    // solving the Ax = f
    // start = high_resolution_clock::now();

    // Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::AMDOrdering<int>> solver;
    // chol_solver.compute(matspec);

    start = high_resolution_clock::now();
    solver.analyzePattern(matspec);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    std::cout << "Time taken to analyze pattern: "
              << duration.count() / 1000000.0 << " seconds" << std::endl;

    start = high_resolution_clock::now();
    solver.factorize(matspec);
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    std::cout << "Time taken to factorize: "
              << duration.count() / 1000000.0 << " seconds" << std::endl;

    // Eigen::IncompleteLUT<double> precond;
    // precond.compute(matspec);
    // Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> solver;
    // Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
    // solver.compute(matspec);
    // solver.preconditioner().setDroptol(std::pow(0.1, 4.0));
    // solver.preconditioner().compute(matspec);

    // solver.compute(matspec);
    start = high_resolution_clock::now();
    vecsol = solver.solve(vecforce);
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    std::cout << "Time taken to set up and solve: "
              << duration.count() / 1000000.0 << " seconds" << std::endl;

    /////////////////////////////////////////////////////////////////////////////////
    start = high_resolution_clock::now();
    // testing class
    auto myinttest = GravityFunctions::Spherical_Integrator<double, EarthModels::PREM>(myprem, vec_numnodes, npoly, 0);
    Eigen::VectorXcd myforce;
    myforce = vecforce;
    Eigen::VectorXcd myclasssol;
    myclasssol = myinttest.solve(myforce);
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    std::cout << "Time taken to set up and solve using my class: "
              << duration.count() / 1000000.0 << " seconds" << std::endl;

    start = high_resolution_clock::now();
    // testing class
    // auto myfullmattest = GravityFunctions::SphericalGeometryPoissonSolver<double, EarthModels::PREM>(myprem, GravityFunctions::CoefficientOrdering::RadialClumped, vec_numnodes, npoly, 20);
    // std::cout << myfullmattest.matrix_length() << std::endl;
    // std::cout << myfullmattest.specmat().nonZeros() << std::endl;
    // Eigen::VectorXcd myforce;
    // myforce = vecforce;
    // Eigen::VectorXcd myclasssol;
    // myclasssol = myinttest.solve(myforce);

    GravityFunctions::PoissonSphericalHarmonic<double, EarthModels::PREM> mybigtest = GravityFunctions::PoissonSphericalHarmonic<double, EarthModels::PREM>(myprem, vec_numnodes, npoly, 0);
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    std::cout << "Time taken to set up and solve using my class: "
              << duration.count() / 1000000.0 << " seconds" << std::endl;

    // finding gravity
    // start = high_resolution_clock::now();

    // find derivative, ie gravity
    Eigen::VectorXd vecderiv = Eigen::VectorXd::Zero(nelem + 1);
    for (int idxelem = 0; idxelem < nelem + 1; ++idxelem)
    {
        if (idxelem < nelem)
        {
            for (int idxpoly = 0; idxpoly < npoly + 1; ++idxpoly)
            {
                vecderiv(idxelem) +=
                    vecsol[idxelem * npoly + idxpoly] * 2.0 / (vec_noderadii[idxelem + 1] - vec_noderadii[idxelem]) * pleg.Derivative(idxpoly, -1.0);
            };
        }
        else
        {
            for (int idxpoly = 0; idxpoly < npoly + 1; ++idxpoly)
            {
                vecderiv(idxelem) +=
                    vecsol[(idxelem - 1) * npoly + idxpoly] * 2.0 / (vec_noderadii[idxelem] - vec_noderadii[idxelem - 1]) * pleg.Derivative(idxpoly, 1.0);
            };
        };
    };
    vecderiv *= 1.0 / rad_scale;

    /////////////////////////////////////////////////////////
    // benchmarking using radial integral
    // start = high_resolution_clock::now();
    // integrals for exact values
    // constant density
    std::vector<double> vec_exactpotential(nelem + 1), vec_exactgravity(nelem + 1);
    for (int idxelem = 0; idxelem < nelem + 1; ++idxelem)
    {
        vec_exactpotential[idxelem] = -2.0 * pi_db * bigg_db * 1000.0 * std::pow(myprem.OuterRadius(), 2) * (1 - 0.333333333333 * std::pow(vec_noderadii[idxelem] * scdiff, 2.0));
        vec_exactgravity[idxelem] = 4.0 * 0.333333333 * pi_db * 1000.0 * bigg_db * myprem.OuterRadius() * vec_noderadii[idxelem] * scdiff;
    }

    ///////////////////////////////////////////////////////////////////

    // PREM exact
    // we need to firstly get the constant part:
    double phi_0;
    phi_0 = 0.0;
    laynum = 0;

    // finding the phi_0 part of the potential
    for (int idxelem = 0; idxelem < nelem; ++idxelem)
    {
        // determining the layer number
        if (!(vec_noderadii[idxelem] < myprem.UpperRadius(laynum) / rad_scale))
        {
            laynum += 1;
        };

        auto funcdens = [&myprem, &rphys, scdiff, laynum, idxelem](double x)
        {
            return rphys(idxelem, x) * myprem.Density(laynum)(rphys(idxelem, x) * scdiff);
            // return rphys(idxelem, x) * 1000.0;
        };
        phi_0 -= 2.0 * pi_db * bigg_db * (vec_noderadii[idxelem + 1] - vec_noderadii[idxelem]) * std::pow(rad_scale, 2.0) * q.Integrate(funcdens);
    }

    // initialising the potential
    std::vector<double> vec_exactprempotential(nelem + 1, 0.0), vec_exactpremgravity(nelem + 1, 0.0);
    vec_exactprempotential[0] += phi_0;

    std::vector<double> vec_massr(nelem + 1, 0.0);
    laynum = 0;
    for (int idxelem = 0; idxelem < nelem; ++idxelem)
    {

        if (!(vec_noderadii[idxelem] < myprem.UpperRadius(laynum) / rad_scale))
        {
            laynum += 1;
        };

        auto funcdens = [&myprem, &rphys, scdiff, laynum, idxelem, &vec_noderadii](double x)
        {
            // return rphys(idxelem, x) * myprem.Density(laynum)(rphys(idxelem, x) * scdiff) * (1 - rphys(idxelem, x) / vec_noderadii[idxelem]);
            return myprem.Density(laynum)(rphys(idxelem, x) * scdiff) * rphys(idxelem, x) * rphys(idxelem, x);
            // return 1000.0 * rphys(idxelem - 1, x) * rphys(idxelem - 1, x);
        };
        vec_massr[idxelem + 1] = vec_massr[idxelem] + 2.0 * pi_db * (vec_noderadii[idxelem + 1] - vec_noderadii[idxelem]) * q.Integrate(funcdens) * std::pow(rad_scale, 3.0);
    };

    std::vector<double> vec_rhorint(nelem + 1, 0.0);
    laynum = 0;
    for (int idxelem = 0; idxelem < nelem; ++idxelem)
    {
        if (!(vec_noderadii[idxelem] < myprem.UpperRadius(laynum) / rad_scale))
        {
            laynum += 1;
        };

        auto funcdens = [&myprem, &rphys, scdiff, laynum, idxelem, &vec_noderadii](double x)
        {
            return rphys(idxelem, x) * myprem.Density(laynum)(rphys(idxelem, x) * scdiff);
            // return 1000.0 * rphys(idxelem - 1, x);
        };
        vec_rhorint[idxelem + 1] = vec_rhorint[idxelem] + 2.0 * pi_db * bigg_db * (vec_noderadii[idxelem + 1] - vec_noderadii[idxelem]) * q.Integrate(funcdens) * std::pow(rad_scale, 2.0);
    };

    // auto stop = high_resolution_clock::now();
    // auto duration = duration_cast<microseconds>(stop - start);
    // std::cout << "Time taken by exact integral solution: "
    //           << duration.count() / 1000000.0 << " seconds" << std::endl;

    ////////////////////////////////////////////////////
    ////////////////////////////////////////////////////
    ////////////////////////////////////////////////////
    ////////////////////////////////////////////////////
    // output to files

    for (int idxelem = 1; idxelem < nelem + 1; ++idxelem)
    {
        vec_exactprempotential[idxelem] += phi_0 - bigg_db * vec_massr[idxelem] / (vec_noderadii[idxelem] * rad_scale) + vec_rhorint[idxelem];
        vec_exactpremgravity[idxelem] = bigg_db * vec_massr[idxelem] / (std::pow(vec_noderadii[idxelem] * rad_scale, 2.0));
    };

    std::string pathtofile = "./work/SpecOutput.out";
    auto file2 = std::ofstream(pathtofile);
    for (int i = 0; i < nelem + 1; ++i)
    {

        file2 << std::setprecision(16) << vec_noderadii[i] << ";"
              << vecsol(i * npoly) << ";" << myclasssol(i * npoly).real() << ";" << vec_exactprempotential[i] << ";" << vecderiv[i] << ";" << vec_exactpremgravity[i] << std::endl;
    };

    // outputting prem:
    int npremx = 1000;
    std::vector<double> vec_premx(npremx, 0.0), vec_premdensity(npremx, 0.0);
    for (int nfill = 0; nfill < npremx; ++nfill)
    {
        vec_premx[nfill] = static_cast<double>(nfill) / static_cast<double>(npremx - 1);
    };
    laynum = 0;
    for (int idxfill = 0; idxfill < npremx; ++idxfill)
    {
        if (!(vec_premx[idxfill] < myprem.UpperRadius(laynum) / rad_scale))
        {
            laynum += 1;
        };
        vec_premdensity[idxfill] = myprem.Density(laynum)(vec_premx[idxfill]);
    };

    laynum = 0;
    pathtofile = "./work/Prem_Density.out";
    auto file = std::ofstream(pathtofile);
    for (int i = 0; i < npremx; ++i)
    {
        if (!(vec_premx[i] < myprem.UpperRadius(laynum) / rad_scale))
        {
            laynum += 1;
        };
        file << std::setprecision(16) << rad_scale / 1000.0 * (1.0 - vec_premx[i]) << ";"
             << vec_premdensity[i] << ";" << myprem.VSV(laynum)(vec_premx[i]) << ";" << myprem.VSH(laynum)(vec_premx[i]) << ";" << myprem.VPV(laynum)(vec_premx[i]) << ";" << myprem.VPH(laynum)(vec_premx[i]) << std::endl;
    };
}
