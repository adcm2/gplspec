#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <GaussQuad/All>
#include <Interpolation/All>
#include <PlanetaryModel/All>
// #include <boost/units/systems/si/current.hpp>
// #include <boost/units/systems/si/electric_potential.hpp>
// #include <boost/units/systems/si/energy.hpp>
// #include <boost/units/systems/si/force.hpp>
// #include <boost/units/systems/si/io.hpp>
// #include <boost/units/systems/si/length.hpp>
// #include <boost/units/systems/si/resistance.hpp>

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
    auto start = high_resolution_clock::now();
    auto myprem = PREM();

    // number of layers and scale factor
    int nlayer = myprem.NumberOfLayers();
    double rad_scale = myprem.OuterRadius();

    // number of nodes within each layer, including those at boundaries
    std::vector<int> vec_numnodes(nlayer, 5);
    vec_numnodes[1] = 5;

    // node values:
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

    // for (int idxtest = 0; idxtest < npoly + 1; ++idxtest)
    // {
    //     std::cout << q.X(idxtest) << std::endl;
    // }

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

    for (int idxelem = 0; idxelem < nelem; ++idxelem)
    {
        for (int idxi = 0; idxi < npoly + 1; ++idxi)
        {
            for (int idxj = 0; idxj < npoly + 1; ++idxj)
            {
                auto funcint = [&func, &func2, &pleg, &idxelem, &idxi, &idxj, &npoly, &rscale, &zetasq](double x)
                {
                    // return rscale(idxelem, x) * rscale(idxelem, x) * func2(idxi, idxj, x) + zetasq * func(idxi, idxj, x);
                    return rscale(idxelem, x) * rscale(idxelem, x) * func2(idxi, idxj, x);
                };
                if (idxi == idxj)
                {
                    tripletList.push_back(T(idxelem * npoly + idxi, idxelem * npoly + idxj, -0.5 * (vec_noderadii[idxelem + 1] - vec_noderadii[idxelem]) * q.Integrate(funcint) + zetasq * q.W(idxi)));
                    // tripletList[idxelem * std::pow((npoly + 1), 2) + idxi * (npoly + 1) + idxj] = T(idxelem * npoly + idxi, idxelem * npoly + idxj, -0.5 * (vec_noderadii[idxelem + 1] - vec_noderadii[idxelem]) * q.Integrate(funcint) + zetasq * q.W(idxi));
                }

                else
                {
                    tripletList.push_back(T(idxelem * npoly + idxi, idxelem * npoly + idxj, -0.5 * (vec_noderadii[idxelem + 1] - vec_noderadii[idxelem]) * q.Integrate(funcint)));
                    // tripletList[idxelem * std::pow((npoly + 1), 2) + idxi * (npoly + 1) + idxj] = T(idxelem * npoly + idxi, idxelem * npoly + idxj, -0.5 * (vec_noderadii[idxelem + 1] - vec_noderadii[idxelem]) * q.Integrate(funcint));
                }
            }
        }
    }

    // final point
    // tripletList[nelem * (npoly + 1) * (npoly + 1) + 1] = T(matlen - 1, matlen - 1, -static_cast<double>(lval + 1) * vec_noderadii[nelem]);
    tripletList.push_back(T(matlen - 1, matlen - 1, -static_cast<double>(lval + 1) * vec_noderadii[nelem]));
    matspec.setFromTriplets(tripletList.begin(), tripletList.end());
    matspec.makeCompressed();

    std::cout << "Size of matrix: " << matlen << std::endl;

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    std::cout << "Setup time pt 1: "
              << duration.count() / 1000000.0 << " seconds" << std::endl;

    start = high_resolution_clock::now();
    double mytest;

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

    // for (int idxi = 0; idxi < npoly + 1; ++idxi)
    // {
    //     Eigen::MatrixXd mat_lderiv = Eigen::MatrixXd::Zero(npoly + 1, npoly + 1);
    //     for (int idxj = 0; idxj < npoly + 1; ++idxj)
    //     {
    //         for (int idxk = 0; idxk < npoly + 1; ++idxk)
    //         {

    //             mat_lderiv(idxk, idxj) = func2(idxi, idxj, q.X(idxk));
    //         }
    //     };
    //     vec_mat_lderiv.push_back(mat_lderiv);
    // };

    // auto functest = [&func, &func2, &pleg, &rscale](int idxelem, int idxi, int idxj, double x)
    // {
    //     return rscale(idxelem, x) * rscale(idxelem, x) * func2(idxi, idxj, x);
    // };
    auto funcnew = [&rphys, &rscale, &vec_mat_lderiv, &q](int idxelem, int idxi, int idxj, int idxk)
    {
        return rphys(idxelem, q.X(idxk)) * rscale(idxelem, q.X(idxk)) * vec_mat_lderiv[idxk](idxi, idxj);
        // return rscale(idxelem, q.X(idxk)) * rscale(idxelem, q.X(idxk)) * vec_mat_lderiv[idxi](idxk, idxj);
    };
    std::vector<double> intvec(npoly + 1);
    for (int idxelem = 0; idxelem < nelem; ++idxelem)
    {

        for (int idxi = 0; idxi < npoly + 1; ++idxi)
        {
            for (int idxj = 0; idxj < npoly + 1; ++idxj)
            {
                // auto funcint = [&func, &func2, &pleg, &idxelem, &idxi, &idxj, &npoly, &rscale, &zetasq](double x)
                // {
                //     // return rscale(idxelem, x) * rscale(idxelem, x) * func2(idxi, idxj, x) + zetasq * func(idxi, idxj, x);
                //     return rscale(idxelem, x) * rscale(idxelem, x) * func2(idxi, idxj, x);
                // };

                for (int idxpoly = 0; idxpoly < npoly + 1; ++idxpoly)
                {
                    // intvec[idxpoly] = functest(idxelem, idxi, idxj, q.X(idxpoly));
                    // intvec[idxpoly] = vec_mat_lderiv[idxpoly](idxi, idxj);
                    intvec[idxpoly] = funcnew(idxelem, idxi, idxj, idxpoly);
                };
                if (idxi == idxj)
                {
                    // mytest = -0.5 * (vec_noderadii[idxelem + 1] - vec_noderadii[idxelem]) * q.Integrate(funcint) + zetasq * q.W(idxi);
                    mytest = std::inner_product(intvec.begin(), intvec.end(), q.Weights().cbegin(), 0.0) + zetasq * q.W(idxi);
                    // tripletList[idxelem * std::pow((npoly + 1), 2) + idxi * (npoly + 1) + idxj] = T(idxelem * npoly + idxi, idxelem * npoly + idxj, -0.5 * (vec_noderadii[idxelem + 1] - vec_noderadii[idxelem]) * q.Integrate(funcint) + zetasq * q.W(idxi));
                }

                else
                {
                    // mytest = -0.5 * (vec_noderadii[idxelem + 1] - vec_noderadii[idxelem]) * q.Integrate(funcint);
                    mytest = std::inner_product(intvec.begin(), intvec.end(), q.Weights().cbegin(), 0.0);
                    // tripletList[idxelem * std::pow((npoly + 1), 2) + idxi * (npoly + 1) + idxj] = T(idxelem * npoly + idxi, idxelem * npoly + idxj, -0.5 * (vec_noderadii[idxelem + 1] - vec_noderadii[idxelem]) * q.Integrate(funcint));
                }
            }
        }
    }

    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    std::cout << "Setup pure integration time: "
              << duration.count() / 1000000.0 << " seconds" << std::endl;

    start = high_resolution_clock::now();
    const double bigg_db = 6.6743 * std::pow(10.0, -11.0);
    const double pi_db = 3.1415926535;
    double scdiff = rad_scale / myprem.OuterRadius();
    double multfact = 2.0 * pi_db * bigg_db * std::pow(rad_scale, 2.0);
    int laynum = 0;
    // std::cout << scdiff << std::endl;
    // force vector
    for (int idxelem = 0; idxelem < nelem; ++idxelem)
    {
        if (!(vec_noderadii[idxelem] < myprem.UpperRadius(laynum) / rad_scale))
        {
            laynum += 1;
        };

        for (int idxpoly = 0; idxpoly < npoly + 1; ++idxpoly)
        {
            auto funcforce = [&myprem, &rphys, &pleg, scdiff, laynum, idxelem, idxpoly](double x)
            {
                return rphys(idxelem, x) * rphys(idxelem, x) * pleg(idxpoly, x) * myprem.Density(laynum)(rphys(idxelem, x) * scdiff);
                // return rphys(idxelem, x) * rphys(idxelem, x) * pleg(idxpoly, x) * 1000.0; // test case with constant density
            };
            vecforce(idxelem * npoly + idxpoly) +=
                multfact * (vec_noderadii[idxelem + 1] - vec_noderadii[idxelem]) * q.Integrate(funcforce);
        };
    };

    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    std::cout << "Setup time pt 2: "
              << duration.count() / 1000000.0 << " seconds" << std::endl;

    ////////////////////////////////////////////////////////////////
    start = high_resolution_clock::now();

    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
    solver.analyzePattern(matspec);
    solver.factorize(matspec);

    // Eigen::IncompleteLUT<double> precond;
    // precond.compute(matspec);
    // Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> solver;
    // Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
    solver.compute(matspec);
    // solver.preconditioner().setDroptol(std::pow(0.1, 4.0));
    // solver.preconditioner().compute(matspec);

    // solver.compute(matspec);

    vecsol = solver.solve(vecforce);
    // vecsol *= 1000.0;

    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    std::cout << "Time taken by solver: "
              << duration.count() / 1000000.0 << " seconds" << std::endl;
    // std::cout << "Number of iterations: " << solver.iterations() << std::endl;

    start = high_resolution_clock::now();
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
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    std::cout << "Time taken by finding derivative: "
              << duration.count() / 1000000.0 << " seconds" << std::endl;
    /////////////////////////////////////////////////////////
    start = high_resolution_clock::now();
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

    for (int idxelem = 0; idxelem < nelem; ++idxelem)
    {
        if (!(vec_noderadii[idxelem] < myprem.UpperRadius(laynum) / rad_scale))
        {
            laynum += 1;
        };
        // std::cout << laynum << std::endl;

        auto funcdens = [&myprem, &rphys, scdiff, laynum, idxelem](double x)
        {
            return rphys(idxelem, x) * myprem.Density(laynum)(rphys(idxelem, x) * scdiff);
            // return rphys(idxelem, x) * 1000.0;
        };
        phi_0 -= 2.0 * pi_db * bigg_db * (vec_noderadii[idxelem + 1] - vec_noderadii[idxelem]) * std::pow(rad_scale, 2.0) * q.Integrate(funcdens);
        // std::cout << phi_0 << std::endl;
    }
    // std::cout << "Exact phi_0 is: " << -2.0 * pi_db * bigg_db * 1000.0 * std::pow(myprem.OuterRadius(), 2.0) << std::endl;
    // std::cout << "phi_0 is: " << phi_0 << std::endl;

    std::vector<double> vec_exactprempotential(nelem + 1, 0.0), vec_exactpremgravity(nelem + 1, 0.0);
    vec_exactprempotential[0] += phi_0;

    std::vector<double> vec_massr(nelem + 1, 0.0);
    laynum = 0;
    for (int idxelem = 1; idxelem < nelem + 1; ++idxelem)
    {
        // if (idxelem - 1 < nelem)
        // {
        if (!(vec_noderadii[idxelem - 1] < myprem.UpperRadius(laynum) / rad_scale))
        {
            laynum += 1;
        };
        // };

        // std::cout << laynum << std::endl;

        auto funcdens = [&myprem, &rphys, scdiff, laynum, idxelem, &vec_noderadii](double x)
        {
            // return rphys(idxelem, x) * myprem.Density(laynum)(rphys(idxelem, x) * scdiff) * (1 - rphys(idxelem, x) / vec_noderadii[idxelem]);
            return myprem.Density(laynum)(rphys(idxelem - 1, x) * scdiff) * rphys(idxelem - 1, x) * rphys(idxelem - 1, x);
            // return 1000.0 * rphys(idxelem - 1, x) * rphys(idxelem - 1, x);
        };
        vec_massr[idxelem] = vec_massr[idxelem - 1] + 2.0 * pi_db * (vec_noderadii[idxelem] - vec_noderadii[idxelem - 1]) * q.Integrate(funcdens) * std::pow(rad_scale, 3.0);
    };

    // std::cout << "Radius is: " << rad_scale << std::endl;
    // std::cout << "Exact mass is: " << 4.0 / 3.0 * pi_db * std::pow(rad_scale, 3.0) * 1000.0 << std::endl;
    // std::cout << "Mass is: " << vec_massr[nelem] << std::endl;

    std::vector<double> vec_rhorint(nelem + 1, 0.0);
    laynum = 0;
    for (int idxelem = 1; idxelem < nelem + 1; ++idxelem)
    {
        // if (idxelem < nelem)
        // {
        if (!(vec_noderadii[idxelem - 1] < myprem.UpperRadius(laynum) / rad_scale))
        {
            laynum += 1;
        };
        // };

        // std::cout << laynum << std::endl;

        auto funcdens = [&myprem, &rphys, scdiff, laynum, idxelem, &vec_noderadii](double x)
        {
            return rphys(idxelem - 1, x) * myprem.Density(laynum)(rphys(idxelem - 1, x) * scdiff);
            // return 1000.0 * rphys(idxelem - 1, x);
        };
        vec_rhorint[idxelem] = vec_rhorint[idxelem - 1] + 2.0 * pi_db * bigg_db * (vec_noderadii[idxelem] - vec_noderadii[idxelem - 1]) * q.Integrate(funcdens) * std::pow(rad_scale, 2.0);
    };

    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    std::cout << "Time taken by exact integral solution: "
              << duration.count() / 1000000.0 << " seconds" << std::endl;

    ////////////////////////////////////////////////////
    // output
    // std::cout << "Exact rhorint: " << 2.0 * pi_db * bigg_db * 1000.0 * std::pow(myprem.OuterRadius(), 2.0) << std::endl;
    // std::cout << "Integrated: " << vec_rhorint[nelem] << std::endl;

    for (int idxelem = 1; idxelem < nelem + 1; ++idxelem)
    {
        vec_exactprempotential[idxelem] += phi_0 - bigg_db * vec_massr[idxelem] / (vec_noderadii[idxelem] * rad_scale) + vec_rhorint[idxelem];
        vec_exactpremgravity[idxelem] = bigg_db * vec_massr[idxelem] / (std::pow(vec_noderadii[idxelem] * rad_scale, 2.0));
    };
    // std::cout << "Hello\n";

    std::string pathtofile = "./work/SpecOutput.out";
    auto file2 = std::ofstream(pathtofile);
    for (int i = 0; i < nelem + 1; ++i)
    {

        file2 << std::setprecision(16) << vec_noderadii[i] << ";"
              << vecsol(i * npoly) << ";" << vec_exactprempotential[i] << ";" << vecderiv[i] << ";" << vec_exactpremgravity[i] << std::endl;
    };
}
