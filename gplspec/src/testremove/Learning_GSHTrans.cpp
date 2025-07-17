#include <FFTWpp/All>
#include <GSHTrans/All>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <concepts>
#include <execution>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <numbers>
#include <random>
#include <TomographyModels/All>
#include <filesystem>
#include <chrono>
#include <PlanetaryModel/All>

using namespace GSHTrans;

using Int = std::ptrdiff_t;

auto RandomDegree()
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<Int> d(3, 7);
    return std::pow(2, d(gen));
}

template <IndexRange NRange>
auto RandomUpperIndex(Int nMax)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    if constexpr (std::same_as<NRange, All>)
    {
        std::uniform_int_distribution<Int> d(-nMax, nMax);
        return d(gen);
    }
    else
    {
        std::uniform_int_distribution<Int> d(0, nMax);
        return d(gen);
    }
}

double azimuthtolongitude(double phi)
{
    return 180.0 / 3.1415926535 * phi;
}
double polartolatitude(double theta)
{
    return 90.0 - (180.0 / 3.1415926535 * theta);
}

int main()
{
    using namespace std::chrono;
    using Real = double;
    using Complex = std::complex<Real>;
    using Scalar = Complex;
    using MRange = All;
    using NRange = NonNegative;
    using Grid = GaussLegendreGrid<Real, MRange, NRange>;
    std::filesystem::path cwd = std::filesystem::current_path();
    std::string fwd = cwd / "modeldata/S40RTS_dvs.nc";
    // std::cout << fwd << std::endl;
    // which *.nc file to use.
    // auto tomo = Tomography("./modeldata/S40RTS_dvs.nc");
    auto tomo = Tomography(fwd);
    // get dvs value at location (depth, lon, lat).

    // finding the deviation in s-wave velocity at depth in mantle with theta and phi
    auto theta = 2.0;
    auto phi = 3.5;
    double thetalat, philon;
    // if (phi > 3.1415926535)
    // {
    //     philon = 180.0 / 3.1415926535 * (-2.0 * 3.1415926535 + phi);
    // }
    // else
    // {
    //     philon = 180.0 / 3.1415926535 * phi;
    // }
    philon = 180.0 / 3.1415926535 * phi;
    thetalat = 90.0 - (180.0 / 3.1415926535 * theta);
    std::cout << tomo.GetValueAt(3000.0, philon, thetalat) << std::endl;
    // std::cout << tomo.GetValueAt(2343.5, 180.0 / 3.1415926535 * phi, thetalat) << std::endl;
    static_assert(PlanetaryModel::SingleParameterDeviationModel<Tomography>);

    // Other APIs.
    std::cout << tomo.GetDepths()[0] << " " << tomo.GetDepths().back() << std::endl;
    std::cout << tomo.GetLongitudes()[0] << " " << tomo.GetLongitudes().back() << std::endl;
    std::cout << tomo.GetLatitudes()[0] << " " << tomo.GetLatitudes().back() << std::endl;
    std::cout << "\n \n " << std::endl;

    // // Let's try find the spherical harmonic decomposition at a particular depth:
    auto start = high_resolution_clock::now();
    auto lMax = 10;
    auto nMax = 0;
    auto grid = Grid(lMax, nMax);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    std::cout << "Time taken to make grid: "
              << duration.count() / 1000000.0 << " seconds" << std::endl;
    // auto n = RandomUpperIndex<NRange>(nMax);

    // Make a random coefficient.
    auto getSize = [](auto lMax, auto n)
    {
        if constexpr (RealFloatingPoint<Scalar>)
        {
            return GSHIndices<NonNegative>(lMax, lMax, n).size();
        }
        else
        {
            return GSHIndices<All>(lMax, lMax, n).size();
        }
    };

    auto size = getSize(lMax, nMax);
    auto flm2 = FFTWpp::vector<Complex>(size, 0.0);
    auto flm = FFTWpp::vector<Complex>(size);
    {
        std::random_device rd{};
        std::mt19937_64 gen{rd()};
        std::normal_distribution<Real> d{0., 1.};
        std::ranges::generate(flm, [&gen, &d]()
                              { return Complex{d(gen), d(gen)}; });

        if constexpr (ComplexFloatingPoint<Scalar>)
        {
            auto flmView = GSHView<Complex, All>(lMax, lMax, nMax, flm.begin());
            flmView(lMax)(lMax) = 0;
        }
        else
        {
            auto flmView = GSHView<Complex, NonNegative>(lMax, lMax, nMax, flm.begin());
            for (auto l : flmView.Degrees())
            {
                flmView(l)(0).imag(0);
            }
            flmView(lMax)(lMax).imag(0);
        }
    }
    flm = flm2;

    flm[3] = 1.0;

    auto f = FFTWpp::vector<Scalar>(grid.NumberOfLongitudes() *
                                    grid.NumberOfCoLatitudes());
    auto glm = FFTWpp::vector<Complex>(size);

    grid.InverseTransformation(lMax, nMax, flm, f);
    grid.ForwardTransformation(lMax, nMax, f, glm);

    // std::cout << std::endl;

    const Complex myi(0.0, 1.0);
    auto fcheck = FFTWpp::vector<Scalar>(grid.NumberOfLongitudes() *
                                         grid.NumberOfCoLatitudes());
    auto ferr = FFTWpp::vector<double>(grid.NumberOfLongitudes() *
                                       grid.NumberOfCoLatitudes());

    for (int idxt = 0; idxt < grid.NumberOfCoLatitudes(); ++idxt)
    {
        for (int idxp = 0; idxp < grid.NumberOfLongitudes(); ++idxp)
        {
            fcheck[idxp + idxt * grid.NumberOfLongitudes()] = -1 / 2.0 * std::sqrt(3.0 / (2.0 * 3.1415926535)) * std::sin(grid.CoLatitudes()[idxt]) * std::exp(myi * grid.LongitudeSpacing() * static_cast<Complex>(idxp));
            // std::cout << idxp + idxt * grid.NumberOfCoLatitudes() << std::endl;
            // std::cout << idxt << " " << idxp << std::endl;
            // fcheck[idxp + idxt * grid.NumberOfLongitudes()] = 1 / 2.0 * std::sqrt(3.0 / (3.1415926535)) * std::cos(grid.CoLatitudes()[idxt]);
        }
    }

    // std::cout << std::endl;
    start = high_resolution_clock::now();
    auto vec_dvs = FFTWpp::vector<Scalar>(grid.NumberOfLongitudes() *
                                          grid.NumberOfCoLatitudes());
    double mydepth = 3000;

    for (int idxt = 0; idxt < grid.NumberOfCoLatitudes(); ++idxt)
    {
        for (int idxp = 0; idxp < grid.NumberOfLongitudes(); ++idxp)
        {
            vec_dvs[idxp + idxt * grid.NumberOfLongitudes()] = tomo.GetValueAt(mydepth, azimuthtolongitude(grid.LongitudeSpacing() * static_cast<double>(idxp)), polartolatitude(grid.CoLatitudes()[idxt]));
        }
    }

    // convert to spherical harmonics:
    auto dvslm = FFTWpp::vector<Complex>(size);
    grid.ForwardTransformation(lMax, nMax, vec_dvs, dvslm);
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    std::cout << "Time taken to perform transforms: "
              << duration.count() / 1000000.0 << " seconds" << std::endl;
    // print:
    // for (const Complex &didx : vec_dvs)
    // {
    //     std::cout << didx << std::endl;
    // }
    // for (const Complex &cidx : dvslm)
    // {
    //     std::cout << cidx << std::endl;
    // }

    for (int myidx = 0; myidx < fcheck.size(); ++myidx)
    {
        ferr[myidx] = std::abs(fcheck[myidx] - f[myidx]);
    }
    auto myerr = std::accumulate(ferr.begin(), ferr.end(), 0.0);

    std::cout << "The error vs exact spherical harmonic is: " << myerr << std::endl;
    // std::cout << "\n";
    std::ranges::transform(flm, glm, flm.begin(),
                           [](auto f, auto g)
                           { return f - g; });

    auto err = *std::ranges::max_element(
        flm, [](auto a, auto b)
        { return std::abs(a) < std::abs(b); });

    std::cout << "Error in forward then inverse is: " << std::abs(err) << std::endl;

    // for (const Complex &cidx : flm)
    // {
    //     std::cout << cidx << " ";
    // };
    // std::cout << "\n";
    // std::cout << flm.size() << std::endl;

    // std::cout << "Colatitudes: \n";
    // auto colatvec = grid.CoLatitudes();
    // for (const Complex &cidx : colatvec)
    // {
    //     std::cout << cidx << " ";
    // };
    // std::cout << std::endl;

    // std::cout << "Longitudes: \n";
    // auto longvec = grid.Longitudes();
    // for (const Complex &cidx : longvec)
    // {
    //     std::cout << cidx << " ";
    // };
    // std::cout << std::endl;
}
