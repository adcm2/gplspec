---
title: "Tutorial 5: Phobos  "
description: "Gravitational potential calculation for Phobos"
layout: default
---

# Phobos
## Outline
This tutorial demonstrates how to construct a model of Phobos. Phobos' density variation is not well constrained. However, its surface topography is relatively well known. The assumption made in the constructor for Phobos can be adapted to a more complex density distribution, but it assumes that the density is homogeneous by default. 

## Implementation
### Information
The model data is taken from Willner 2014. The surface radius is given in terms of a real spherical harmonic decomposition. The constructor used here assumes that the model data is in this format. The constructor needs the density, a path to the file and some more parameters. 

### Rotation
We also introduce the rotational functionality. The way that this is implemented is one specifies two points which are rotated onto the equator, with the first situated at the prime meridian. The potential can be outputted either in the original orientation or in a rotated manner. 

### Source code
```cpp
#include <gplspec/All>

int
main() {
   using namespace GeneralEarthModels;
   using namespace Gravity_Tools;
   using namespace SimpleModels;
   using FLOAT = double;
   std::string pathtofile = "modeldata/TABLEA1.DAT";

   // parameters
   double physdensity = 1860.0, lengthnorm = 10000.0,  usedensity = 1342.0,
    massnorm = std::pow(lengthnorm, 3.0) * usedensity, timenorm = 3600.0,
     scaledensity = physdensity / massnorm * std::pow(lengthnorm, 3.0);
   int lMax = 100;
   int npoly = 5;

   // constructing model
   Density3D phobos(scaledensity, pathtofile, npoly, lMax, lengthnorm, timenorm,
                    massnorm, 0.1, 1.5);

   // get gravitational field
   auto stdvec_potsol =
       FindGravitationalPotential(phobos, std::pow(10.0, -12.0));

   auto stdvec_potsol2 =
       FindGravitationalPotential(phobos, std::pow(10.0, -3.0));
       
 
   // coordinates for rotation
   double theta1, theta2;
   double phi1, phi2;

   // vectors containing theta,phi information
   std::vector<double> vec_ang1{theta1, phi1}, vec_ang2{theta2, phi2};

   //////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////

   // output test
   std::string pathtofolder1 = "./work/Phobos/Cartesian";
   std::string pathtofolder2 = "./work/Phobos/Spherical";
   std::string pathtofolder3 = "./work/Phobos/Sensitivity";
   std::string pathtofile1 =
       pathtofolder2 + "/MatrixSolutionReferentialRotated.out";
   std::string pathtofile2 = pathtofolder2 + "/MatrixSolutionPhysical.out";
   std::string pathtofile3 = pathtofolder2 + "/MatrixSolutionRotated.out";
   std::string pathtofile4 =
       pathtofolder2 + "/PhysicalDensitySolutionRotated.out";
   std::string pathtofile5 =
       pathtofolder2 + "/ReferentialDensitySolutionRotated.out";

   phobos.CartesianOutputAtElement(pathtofolder1, stdvec_potsol);
   phobos.PhysicalOutputAtElement(pathtofolder2, stdvec_potsol);
   phobos.ReferentialOutputAtElement(pathtofolder3, stdvec_senskernel, true);
   phobos.ReferentialOutputSlice(pathtofolder3, stdvec_senskernel, true);
   phobos.PhysicalOutputSlice(pathtofolder2, stdvec_potsol);
   phobos.ReferentialOutputRotated(pathtofile1, vec_ang1, vec_ang2,
                                   stdvec_potsol);
   phobos.PhysicalOutputRotated(pathtofile2, vec_ang12, vec_ang22,
                                stdvec_potsol);

   return 0;
}


```





