#include <Eigen/Core>
#include <Eigen/Dense>
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

    auto myprem = PREM();

    // number of layers and scale factor
    int nlayer = myprem.NumberOfLayers();
    double rad_scale = myprem.OuterRadius();

    // polynomial order
    int npoly = 6;
    std::cout << "Order of polynomial: \n";
    std::cin >> npoly;
    auto start = high_resolution_clock::now();
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

    // int npoly = 8;                  // number of polynomials within each element
    int matlen = nelem * npoly + 1; // size of matrix

    // function:
    int N = 5;
    std::vector<double> x(N), y(N);

    // // finite element matrices
    Eigen::MatrixXd matspec = Eigen::MatrixXd::Zero(matlen, matlen);
    Eigen::VectorXd vecforce = Eigen::VectorXd::Zero(matlen);
    Eigen::VectorXd vecsol = Eigen::VectorXd::Zero(matlen);

    // generate Gauss grid
    auto q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(npoly + 1);
    // std::vector<double> vval(nelem, 1.0);

    // // vector of Lagrange polynomials:

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
    for (int idxelem = 0; idxelem < nelem; ++idxelem)
    {
        for (int idxi = 0; idxi < npoly + 1; ++idxi)
        {
            for (int idxj = 0; idxj < npoly + 1; ++idxj)
            {
                auto funcint = [&func, &func2, pleg, idxelem, idxi, idxj, npoly, &rscale, lval](double x)
                {
                    return rscale(idxelem, x) * rscale(idxelem, x) * func2(idxi, idxj, x) + static_cast<double>(lval) * static_cast<double>(lval + 1) * func(idxi, idxj, x);
                };
                matspec(idxelem * npoly + idxi, idxelem * npoly + idxj) -= 0.5 * (vec_noderadii[idxelem + 1] - vec_noderadii[idxelem]) * q.Integrate(funcint);
            }
        }
    }

    // final point
    matspec(matlen - 1, matlen - 1) -= static_cast<double>(lval + 1) * vec_noderadii[nelem];

    std::cout << "Size of matrix: " << matlen << std::endl;

    const double bigg_db = 6.6743 * std::pow(10.0, -11.0);
    const double pi_db = 3.1415926535;
    double scdiff = rad_scale / myprem.OuterRadius();
    double multfact = 2.0 * pi_db * bigg_db * std::pow(rad_scale, 2.0);
    int laynum = 0;
    std::cout << scdiff << std::endl;
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
    }

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    std::cout << "Time taken for setup: "
              << duration.count() / 1000000.0 << " seconds" << std::endl;

    start = high_resolution_clock::now();
    Eigen::FullPivLU<Eigen::MatrixXd> solver(matspec);
    vecsol = solver.solve(vecforce);
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    std::cout << "Time taken by solver: "
              << duration.count() / 1000000.0 << " seconds" << std::endl;
    // vecsol *= 1000.0;

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

    std::cout << "Radius is: " << rad_scale << std::endl;
    std::cout << "Exact mass is: " << 4.0 / 3.0 * pi_db * std::pow(rad_scale, 3.0) * 1000.0 << std::endl;
    std::cout << "Mass is: " << vec_massr[nelem] << std::endl;

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
