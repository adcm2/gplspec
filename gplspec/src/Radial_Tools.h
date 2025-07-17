#ifndef RADIAL_TOOLS_H
#define RADIAL_TOOLS_H

#include <GaussQuad/All>
#include <PlanetaryModel/All>

namespace Radial_Tools {

// returns mapping of x from [-1,1] to [x1,x2]
template <typename FLOAT>
FLOAT
StandardIntervalMap(const FLOAT &x, const FLOAT &x1, const FLOAT &x2) {
   return ((x2 - x1) * x + (x1 + x2)) * 0.5;
};

// determines radii of all nodes within the decomposition
template <typename FLOAT>
std::vector<FLOAT>
All_Node(const std::vector<FLOAT> &vec_noderadii,
         const GaussQuad::Quadrature1D<FLOAT> &q) {

   // IMPORTANT: MAX ELEMENT SIZE MUST BE IN SAME UNITS AS MODEL//

   // sizes of various vectors:
   std::size_t nelem = vec_noderadii.size() - 1;
   int npoly = q.N() - 1;
   // std::cout << "npoly: " << npoly << "\n";
   std::size_t matlen = nelem * npoly + 1;    // size of vector
   std::vector<FLOAT> vec_allradii(matlen);   // declare vector
   for (int idx = 0; idx < nelem; ++idx) {
      for (int idx2 = 0; idx2 < npoly; ++idx2) {
         vec_allradii[idx * npoly + idx2] = StandardIntervalMap(
             q.X(idx2), vec_noderadii[idx], vec_noderadii[idx + 1]);
      }
   }
   vec_allradii[0] = 0.0;
   vec_allradii[matlen - 1] = vec_noderadii[nelem];
   return vec_allradii;
};

// determines radii of all nodes within the decomposition
template <typename FLOAT>
std::vector<std::vector<FLOAT>>
All_Node_New(const std::vector<FLOAT> &vec_noderadii,
             const GaussQuad::Quadrature1D<FLOAT> &q) {

   // IMPORTANT: MAX ELEMENT SIZE MUST BE IN SAME UNITS AS MODEL//

   // sizes of various vectors:
   std::size_t nelem = vec_noderadii.size() - 1;
   int npoly = q.N() - 1;
   // std::cout << "npoly: " << npoly << "\n";
   std::vector<std::vector<FLOAT>> vec_allradii(nelem);   // declare vector

   // fill out
   for (int idx = 0; idx < nelem; ++idx) {
      auto tmp_radius = std::vector<FLOAT>(npoly + 1, 0.0);
      for (int idx2 = 0; idx2 < npoly + 1; ++idx2) {
         tmp_radius[idx2] = StandardIntervalMap(q.X(idx2), vec_noderadii[idx],
                                                vec_noderadii[idx + 1]);
      }
      vec_allradii[idx] = tmp_radius;
   }
   // vec_allradii[0] = 0.0;
   // vec_allradii[matlen - 1] = vec_noderadii[nelem];
   return vec_allradii;
};

// class which contains nodes as well as info on indices where discontinuity
// occurs

class node_info {
 public:
   node_info() {};   // default

   template <class sphericalmodel>
      requires PlanetaryModel::BasicSphericalDensityModel<sphericalmodel, int,
                                                          double>
   node_info(const sphericalmodel &, const GaussQuad::Quadrature1D<double> &,
             const double, const double);

   // outputs
   auto idx_discont() { return _vec_discontinuity_indices; };
   auto ElementNodes() { return _vec_noderadii; };
   auto AllNodes() { return _vec_allnodes; };
   auto NumberOfNodes() { return _vec_noderadii.size(); };
   auto NumberOfElements() { return _vec_noderadii.size() - 1; };
   auto TotalNumberOfNodes() { return _vec_allnodes.size(); };
   auto Node(int idx) {
      assert(idx < _vec_noderadii.size() && "Not in range");
      return _vec_noderadii[idx];
   };
   auto AllNode(int idx) const {
      assert(idx < _vec_allnodes.size() && "Not in range");
      return _vec_allnodes[idx];
   };
   auto NodeRadius(int idxelem, int idxpoly) {
      return _vvec_allnodes[idxelem][idxpoly];
   };
   auto LayerNumber(int idxelem) { return _vec_layers[idxelem]; };

   auto IdxLowerRadius(int idx) { return _vec_discontinuity_indices[idx]; };
   auto IdxUpperRadius(int idx) { return _vec_discontinuity_indices[idx + 1]; };
   auto NumberOfLayers() const {
      return _vec_discontinuity_indices.size() - 1;
   };
   auto LowerRadius(int idx) {
      return _vec_noderadii[_vec_discontinuity_indices[idx]];
   };
   auto UpperRadius(int idx) {
      return _vec_noderadii[_vec_discontinuity_indices[idx + 1]];
   };
   auto OuterRadius() { return _vec_noderadii.back(); };

 private:
   std::vector<int> _vec_discontinuity_indices, _vec_layers;
   std::vector<double> _vec_noderadii, _vec_allnodes;
   std::vector<std::vector<double>> _vvec_allnodes;
};

// setting
template <class sphericalmodel>
   requires PlanetaryModel::BasicSphericalDensityModel<sphericalmodel, int,
                                                       double>
node_info::node_info(const sphericalmodel &mymodel,
                     const GaussQuad::Quadrature1D<double> &q,
                     const double max_element_size, const double ballrad) {

   // first element of both
   _vec_noderadii.push_back(0.0);
   _vec_discontinuity_indices.push_back(0);

   for (int idx = 0; idx < mymodel.NumberOfLayers(); ++idx) {
      // Determining the spacing of nodes within each layer
      double laydepth = mymodel.UpperRadius(idx) - mymodel.LowerRadius(idx);
      int numelements = std::ceil(laydepth / max_element_size);

      // push back to index vector
      _vec_discontinuity_indices.push_back(_vec_discontinuity_indices.back() +
                                           numelements);
      double nodespacing = laydepth / static_cast<double>(numelements);

      // filling out the node radii vector
      for (int idxelem = 0; idxelem < numelements; ++idxelem) {
         if (idxelem < numelements - 1) {
            double currentval = _vec_noderadii.back();
            _vec_noderadii.push_back(currentval + nodespacing);
            _vec_layers.push_back(idx);
         } else {
            _vec_noderadii.push_back(mymodel.UpperRadius(idx));
            _vec_layers.push_back(idx);
         }
      }
   };

   // adding in the outer elements between the edge of the planet and the
   {
      // get size and number of elements
      double laydepth = ballrad - mymodel.OuterRadius();
      assert(laydepth > 0.0 &&
             "Ball radius must be larger than outer radius of planet");
      int numelem = std::ceil(laydepth / max_element_size);
      double nodespacing = laydepth / static_cast<double>(numelem);

      // finish filling out vector
      for (int idxelem = 0; idxelem < numelem; ++idxelem) {
         if (idxelem < numelem - 1) {
            double currentval = _vec_noderadii.back();
            _vec_noderadii.push_back(currentval + nodespacing);
            _vec_layers.push_back(mymodel.NumberOfLayers());
         } else {
            _vec_noderadii.push_back(ballrad);
            _vec_layers.push_back(mymodel.NumberOfLayers());
         }
      }
   }

   // all nodes
   _vec_allnodes = All_Node(_vec_noderadii, q);
   _vvec_allnodes = All_Node_New(_vec_noderadii, q);
   // check vvec
   // std::cout << "length: " << _vvec_allnodes.size()
   //           << ", of individual: " << _vvec_allnodes[0].size() << "\n";
   // std::cout << "nelem: " << _vec_noderadii.size() - 1 << ", npoly: " <<
   // q.N()
   //           << "\n";
};

// class which contains nodes as well as info on indices where discontinuity
// occurs

class RadialMesh {
 public:
   // constructors
   RadialMesh() {};   // default

   // simple constructor
   RadialMesh(const double, const double, const double,
              const GaussQuad::Quadrature1D<double> &);

   template <class sphericalmodel>
      requires PlanetaryModel::BasicSphericalDensityModel<sphericalmodel, int,
                                                          double>
   RadialMesh(const sphericalmodel &, const GaussQuad::Quadrature1D<double> &,
              const double, const double, bool = true);

   // return functions of nodes
   auto NumberOfElements() const { return _vec_nodes.size(); };
   auto NumNodesWithinElement() const { return _vec_nodes[0].size(); };
   auto NodeRadius(int idxelem, int idxnode) const {
      return _vec_nodes[idxelem][idxnode];
   };

   // information on elements
   auto LayerNumber(int idxelem) const { return _vec_layers[idxelem]; };
   auto NumberOfLayers() const { return _vec_layers.back() + 1; };
   auto ElementLowerRadius(int idxelem) const {
      return _vec_nodes[idxelem][0];
   };
   auto ElementUpperRadius(int idxelem) const {
      return _vec_nodes[idxelem].back();
   };
   auto ElementWidth(int idxelem) const {
      return _vec_nodes[idxelem].back() - _vec_nodes[idxelem][0];
   };
   auto OuterRadius() const { return _vec_nodes.back().back(); };
   auto PlanetRadius() const { return planet_radius; };

 private:
   std::vector<int> _vec_layers;
   std::vector<std::vector<double>> _vec_nodes;
   double planet_radius;
};

RadialMesh::RadialMesh(const double planetradius, const double ballrad,
                       const double max_element_size,
                       const GaussQuad::Quadrature1D<double> &q)
    : planet_radius{planetradius} {
   // lambda to generate nodes between two radii
   auto nodegenerate = [&q](double lowerradius, double upperradius) {
      // vector to contain all nodes in element
      auto tmp_radius = std::vector<double>(q.N(), 0.0);
      tmp_radius[0] = lowerradius;
      for (int idx2 = 1; idx2 < q.N() - 1; ++idx2) {
         tmp_radius[idx2] =
             StandardIntervalMap(q.X(idx2), lowerradius, upperradius);
      }
      tmp_radius[q.N() - 1] = upperradius;
      return tmp_radius;
   };
   std::vector<double> element_width;
   std::vector<int> num_elements_layer;
   num_elements_layer.push_back(std::ceil(planet_radius / max_element_size));
   element_width.push_back(planet_radius /
                           static_cast<double>(num_elements_layer[0]));
   num_elements_layer.push_back(
       std::ceil((ballrad - planet_radius) / max_element_size));
   element_width.push_back((ballrad - planet_radius) /
                           static_cast<double>(num_elements_layer[1]));

   // number of elements in total
   auto totlength =
       std::reduce(num_elements_layer.begin(), num_elements_layer.end());
   _vec_nodes.resize(totlength);
   _vec_layers.resize(totlength);

   //////////////////////////////////////////////////////////////////////
   // fill out _vec_nodes
   int idx_global = 0;
   // filling up to outer edge of planet
   // for (int idx = 0; idx < totlayers; ++idx) {
   // up to final element
   for (int idxelem = 0; idxelem < num_elements_layer[0] - 1; ++idxelem) {
      auto lowerradius = idxelem * element_width[0];
      auto upperradius = lowerradius + element_width[0];
      _vec_nodes[idx_global] = nodegenerate(lowerradius, upperradius);
      _vec_layers[idx_global] = 0;
      ++idx_global;
   }

   // getting correct upper radius
   {
      auto lowerradius = planetradius - element_width[0];
      auto upperradius = planetradius;
      _vec_nodes[idx_global] = nodegenerate(lowerradius, upperradius);
      _vec_layers[idx_global] = 0;
      ++idx_global;
   }
   // };

   // filling layer from outer edge to ball
   {

      // first part
      for (int idxelem = 0; idxelem < num_elements_layer[1] - 1; ++idxelem) {
         auto lowerradius = planetradius + idxelem * element_width[1];
         auto upperradius = lowerradius + element_width[1];
         _vec_nodes[idx_global] = nodegenerate(lowerradius, upperradius);
         _vec_layers[idx_global] = 1;
         ++idx_global;
      }

      // second part
      {
         auto lowerradius = ballrad - element_width[1];
         auto upperradius = ballrad;
         _vec_nodes[idx_global] = nodegenerate(lowerradius, upperradius);
         _vec_layers[idx_global] = 1;
         ++idx_global;
      }
   }
};

// RadialMesh::RadialMesh(const double planetradius, const double ballrad,
//                        const double max_element_size,
//                        const GaussQuad::Quadrature1D<double> &q)
//     : planet_radius{planetradius} {
//    // lambda to generate nodes between two radii
//    auto nodegenerate = [&q](double lowerradius, double upperradius) {
//       // vector to contain all nodes in element
//       auto tmp_radius = std::vector<double>(q.N(), 0.0);
//       tmp_radius[0] = lowerradius;
//       for (int idx2 = 1; idx2 < q.N() - 1; ++idx2) {
//          tmp_radius[idx2] =
//              StandardIntervalMap(q.X(idx2), lowerradius, upperradius);
//       }
//       tmp_radius[q.N() - 1] = upperradius;
//       return tmp_radius;
//    };
//    std::vector<double> element_width;
//    std::vector<int> num_elements_layer;
//    // if(max_element_size > 0.01){
//    num_elements_layer.push_back(1);
//    element_width.push_back(planet_radius * 0.01);
//    // }
//    num_elements_layer.push_back(
//        std::ceil(planet_radius * 0.99 / max_element_size));
//    element_width.push_back(planet_radius * 0.99 /
//                            static_cast<double>(num_elements_layer[1]));
//    num_elements_layer.push_back(
//        std::ceil((ballrad - planet_radius) / max_element_size));
//    element_width.push_back((ballrad - planet_radius) /
//                            static_cast<double>(num_elements_layer[2]));
//    std::vector<double> rad_elements;
//    rad_elements.push_back(0.0);
//    rad_elements.push_back(planet_radius * 0.01);
//    rad_elements.push_back(planet_radius);
//    rad_elements.push_back(ballrad);

//    // std::cout << "\n";
//    // for (auto idx : element_width) {
//    //    std::cout << idx << "\n";
//    // }
//    // std::cout << "\n";
//    // std::cout << "\n";
//    // for (auto idx : num_elements_layer) {
//    //    std::cout << idx << "\n";
//    // }
//    // std::cout << "\n";
//    // number of elements in total
//    auto totlength =
//        std::reduce(num_elements_layer.begin(), num_elements_layer.end());
//    _vec_nodes.resize(totlength);
//    _vec_layers.resize(totlength);
//    // std::cout << "Hello\n";
//    //////////////////////////////////////////////////////////////////////
//    // fill out _vec_nodes
//    int idx_global = 0;
//    // filling up to outer edge of planet
//    // for (int idx = 0; idx < totlayers; ++idx) {
//    // up to final element
//    for (int idxtest = 0; idxtest < 2; ++idxtest) {
//       for (int idxelem = 0; idxelem < num_elements_layer[idxtest] - 1;
//            ++idxelem) {
//          auto lowerradius =
//              rad_elements[idxtest] + idxelem * element_width[idxtest];
//          auto upperradius = lowerradius + element_width[idxtest];
//          _vec_nodes[idx_global] = nodegenerate(lowerradius, upperradius);
//          _vec_layers[idx_global] = idxtest;
//          ++idx_global;
//          // if (idxelem == 0) {
//          //    std::cout << "\n";
//          //    std::cout << lowerradius << " " << upperradius << "\n\n";
//          // }
//       }

//       {
//          double lowerradius;
//          if (num_elements_layer[idxtest] == 1) {
//             lowerradius = rad_elements[idxtest];
//          } else {
//             lowerradius = rad_elements[idxtest + 1] - element_width[idxtest];
//          }
//          auto upperradius = rad_elements[idxtest + 1];
//          _vec_nodes[idx_global] = nodegenerate(lowerradius, upperradius);
//          _vec_layers[idx_global] = idxtest;
//          ++idx_global;
//       }
//    }

//    // getting correct upper radius

//    // };

//    // filling layer from outer edge to ball
//    {

//       // first part
//       for (int idxelem = 0; idxelem < num_elements_layer[2] - 1; ++idxelem) {
//          auto lowerradius = planetradius + idxelem * element_width[2];
//          auto upperradius = lowerradius + element_width[2];
//          _vec_nodes[idx_global] = nodegenerate(lowerradius, upperradius);
//          _vec_layers[idx_global] = 2;
//          ++idx_global;
//       }

//       // second part
//       {
//          auto lowerradius = ballrad - element_width[2];
//          auto upperradius = ballrad;
//          _vec_nodes[idx_global] = nodegenerate(lowerradius, upperradius);
//          _vec_layers[idx_global] = 2;
//          ++idx_global;
//       }
//    }
//    // std::cout << "Hello 2\n";
//    // for (auto idx : _vec_layers) {
//    //    std::cout << idx << "\n";
//    // }
//    // std::cout << "\n";
//    // for (auto idx : _vec_nodes) {
//    //    for (auto idx2 : idx) {
//    //       std::cout << idx2 << "\n";
//    //    }
//    // }
//    // std::cout << "\n";
//    // std::cout << "\n";
// };

// setting
template <class sphericalmodel>
   requires PlanetaryModel::BasicSphericalDensityModel<sphericalmodel, int,
                                                       double>
RadialMesh::RadialMesh(const sphericalmodel &mymodel,
                       const GaussQuad::Quadrature1D<double> &q,
                       const double max_element_size, const double ballrad,
                       bool incball)
    : planet_radius(mymodel.OuterRadius()) {

   // lambda to generate nodes between two radii
   auto nodegenerate = [&q](double lowerradius, double upperradius) {
      // vector to contain all nodes in element
      auto tmp_radius = std::vector<double>(q.N(), 0.0);
      tmp_radius[0] = lowerradius;
      for (int idx2 = 1; idx2 < q.N() - 1; ++idx2) {
         tmp_radius[idx2] =
             StandardIntervalMap(q.X(idx2), lowerradius, upperradius);
      }
      tmp_radius[q.N() - 1] = upperradius;
      return tmp_radius;
   };

   // total layers and vectors to store data on elements
   int totlayers = mymodel.NumberOfLayers();
   int numlayers;
   if (incball) {
      numlayers = totlayers + 1;
   } else {
      numlayers = totlayers;
   }
   std::vector<double> element_width(numlayers);
   std::vector<int> num_elements_layer(numlayers);

   // finding width of elements and number of elements within each layer
   for (int idx = 0; idx < totlayers; ++idx) {
      double laydepth = mymodel.UpperRadius(idx) - mymodel.LowerRadius(idx);
      num_elements_layer[idx] = std::ceil(laydepth / max_element_size);
      element_width[idx] =
          laydepth / static_cast<double>(num_elements_layer[idx]);
   };

   if (incball) {
      double laydepth = ballrad - mymodel.OuterRadius();
      num_elements_layer[totlayers] = std::ceil(laydepth / max_element_size);
      element_width[totlayers] =
          laydepth / static_cast<double>(num_elements_layer[totlayers]);
   }

   // number of elements in total
   auto totlength =
       std::reduce(num_elements_layer.begin(), num_elements_layer.end());
   _vec_nodes.resize(totlength);
   _vec_layers.resize(totlength);

   //////////////////////////////////////////////////////////////////////
   // fill out _vec_nodes
   int idx_global = 0;
   // filling up to outer edge of planet
   for (int idx = 0; idx < totlayers; ++idx) {
      // up to final element
      for (int idxelem = 0; idxelem < num_elements_layer[idx] - 1; ++idxelem) {
         auto lowerradius =
             mymodel.LowerRadius(idx) + idxelem * element_width[idx];
         auto upperradius = lowerradius + element_width[idx];
         _vec_nodes[idx_global] = nodegenerate(lowerradius, upperradius);
         _vec_layers[idx_global] = idx;
         ++idx_global;
      }

      // getting corret upper radius
      {
         auto lowerradius = mymodel.UpperRadius(idx) - element_width[idx];
         auto upperradius = mymodel.UpperRadius(idx);
         _vec_nodes[idx_global] = nodegenerate(lowerradius, upperradius);
         _vec_layers[idx_global] = idx;
         ++idx_global;
      }
   };

   // filling layer from outer edge to ball
   if (incball) {

      // first part
      for (int idxelem = 0; idxelem < num_elements_layer[totlayers] - 1;
           ++idxelem) {
         auto lowerradius =
             mymodel.OuterRadius() + idxelem * element_width[totlayers];
         auto upperradius = lowerradius + element_width[totlayers];
         _vec_nodes[idx_global] = nodegenerate(lowerradius, upperradius);
         _vec_layers[idx_global] = totlayers;
         ++idx_global;
      }

      // second part
      {
         auto lowerradius = ballrad - element_width[totlayers];
         auto upperradius = ballrad;
         _vec_nodes[idx_global] = nodegenerate(lowerradius, upperradius);
         _vec_layers[idx_global] = totlayers;
         ++idx_global;
      }
   }
};

// determines radii of edges of elements, with a maximum element size specified
template <typename FLOAT, template <typename, typename> class sphericalmodel>
   requires PlanetaryModel::SphericalElasticModel<sphericalmodel<FLOAT, int>>
std::vector<FLOAT>
Radial_Node(const sphericalmodel<FLOAT, int> &mymodel,
            const FLOAT max_element_size) {

   // IMPORTANT: MAX ELEMENT SIZE MUST BE IN SAME UNITS AS MODEL//

   // declare vec_noderadii
   std::vector<FLOAT> vec_noderadii;   // node radii

   vec_noderadii.push_back(0.0);

   // we want to have normalised radii. Thus define normalisation:
   // FLOAT mynorm = 1.0 / mymodel.OuterRadius();
   // FLOAT mynorm = 1.0 / mymodel.LengthNorm();
   // FLOAT normmaxsize = mynorm * max_element_size;

   for (int idx = 0; idx < mymodel.NumberOfLayers(); ++idx) {
      // Determining the spacing of nodes within each layer
      FLOAT laydepth = mymodel.UpperRadius(idx) - mymodel.LowerRadius(idx);
      int numelements = std::ceil(laydepth / max_element_size);
      FLOAT nodespacing = laydepth / static_cast<FLOAT>(numelements);

      // filling out the node radii vector
      for (int idxelem = 0; idxelem < numelements; ++idxelem) {
         if (idxelem < numelements - 1) {
            FLOAT currentval = vec_noderadii.back();
            vec_noderadii.push_back(currentval + nodespacing);
         } else {
            vec_noderadii.push_back(mymodel.UpperRadius(idx));
         }
      }
   };
   return vec_noderadii;
};

// determines radii of edges of elements, with a maximum element size specified
template <typename FLOAT, template <typename, typename> class sphericalmodel>
   requires PlanetaryModel::SphericalElasticModel<sphericalmodel<FLOAT, int>>
std::vector<FLOAT>
Radial_Node(const sphericalmodel<FLOAT, int> &mymodel,
            const FLOAT max_element_size, const FLOAT ballrad) {

   // IMPORTANT: MAX ELEMENT SIZE MUST BE IN SAME UNITS AS MODEL//

   // declare vec_noderadii
   std::vector<FLOAT> vec_noderadii;   // node radii

   vec_noderadii.push_back(0.0);

   // we want to have normalised radii. Thus define normalisation:
   // FLOAT mynorm = 1.0 / mymodel.LengthNorm();
   FLOAT normmaxsize = max_element_size;

   for (int idx = 0; idx < mymodel.NumberOfLayers(); ++idx) {
      // Determining the spacing of nodes within each layer
      FLOAT laydepth = mymodel.UpperRadius(idx) - mymodel.LowerRadius(idx);
      int numelements = std::ceil(laydepth / normmaxsize);
      FLOAT nodespacing = laydepth / static_cast<FLOAT>(numelements);

      // filling out the node radii vector
      for (int idxelem = 0; idxelem < numelements; ++idxelem) {
         if (idxelem < numelements - 1) {
            FLOAT currentval = vec_noderadii.back();
            vec_noderadii.push_back(currentval + nodespacing);
         } else {
            vec_noderadii.push_back(mymodel.UpperRadius(idx));
         }
      }
   };

   // adding in the outer elements between the edge of the planet and the ball
   // FLOAT multfact = 1.5;
   FLOAT ballrad2 = ballrad;
   // std::cout << "ball radius: " << ballrad2 << std::endl;
   {
      // get size and number of elements
      FLOAT laydepth = ballrad2 - mymodel.OuterRadius();
      int numelem = std::ceil(laydepth / normmaxsize);
      FLOAT nodespacing = laydepth / static_cast<FLOAT>(numelem);

      // finish filling out vector
      for (int idxelem = 0; idxelem < numelem; ++idxelem) {
         if (idxelem < numelem - 1) {
            FLOAT currentval = vec_noderadii.back();
            vec_noderadii.push_back(currentval + nodespacing);
         } else {
            vec_noderadii.push_back(ballrad2);
         }
      }
   }
   return vec_noderadii;
};

}   // namespace Radial_Tools

// first element of both
// _vec_noderadii.push_back(0.0);
// _vec_discontinuity_indices.push_back(0);
// int totelem = 0;
// bool NotFirstElement{false};

// auto newnodes = [&q](double lowerradius, double upperradius,
//                      double maxstep) {
//    double laydepth = upperradius - lowerradius;
//    int numelements = std::ceil(laydepth / maxstep);

//    // element spacing
//    double elementwidth = laydepth / static_cast<double>(numelements);

//    // bool NotFirstElement{false};
//    std::vector<std::vector<double>> _vec_nodes;

//    // filling out the node radii vector
//    for (int idxelem = 0; idxelem < numelements; ++idxelem) {
//       // lower radius
//       double elemlowerradius =
//           lowerradius + static_cast<double>(idxelem) * elementwidth;

//       // upper radius
//       double elemupperradius;
//       if (idxelem < numelements - 1) {
//          elemupperradius = elemlowerradius + elementwidth;
//       } else {
//          elemupperradius = upperradius;
//       }

//       // vector to contain all nodes in element
//       auto tmp_radius = std::vector<double>(q.N(), 0.0);
//       tmp_radius[0] = elemlowerradius;
//       for (int idx2 = 1; idx2 < q.N() - 1; ++idx2) {
//          tmp_radius[idx2] = StandardIntervalMap(q.X(idx2),
//          elemlowerradius,
//                                                 elemupperradius);
//       }
//       tmp_radius[q.N() - 1] = elemupperradius;

//       // nodes within element
//       _vec_nodes.push_back(tmp_radius);

//       // if (!NotFirstElement) {
//       //    NotFirstElement = true;
//       // }
//    }
//    return _vec_nodes;
// };

//////////////////////////////////////////
// for (int idx = 0; idx < mymodel.NumberOfLayers(); ++idx) {
//    // elements within layer
//    auto tmp_nodes = newnodes(mymodel.LowerRadius(idx),
//                              mymodel.UpperRadius(idx), max_element_size);
//    std::copy(tmp_nodes.begin(), tmp_nodes.end(),
//              _vec_nodes.begin() + layernumelements[idx]);
// }

// layer between edge of planet and ball
// auto tmp_nodes = newnodes(mymodel.OuterRadius(), ballrad,
// max_element_size); std::copy(tmp_nodes.begin(), tmp_nodes.end(),
//           _vec_nodes.begin() + layernumelements.back());

////////////////////////////////////////////

// layer number for element
// int numelements = tmp_nodes.size();
// std::vector<int> tmp_layers(numelements, mymodel.NumberOfLayers());
// _vec_layers.insert(_vec_layers.end(), tmp_layers.begin(),
// tmp_layers.end());

// for (int idx = 0; idx < mymodel.NumberOfLayers(); ++idx) {
//    // Determining the spacing of nodes within each layer
//    double laydepth = mymodel.UpperRadius(idx) - mymodel.LowerRadius(idx);
//    int numelements = std::ceil(laydepth / max_element_size);

//    // element spacing
//    double elementwidth = laydepth / static_cast<double>(numelements);

//    // filling out the node radii vector
//    for (int idxelem = 0; idxelem < numelements; ++idxelem) {
//       // lower radius
//       double lowerradius;
//       if (NotFirstElement) {
//          lowerradius = _vec_nodes.back().back();
//       } else {
//          lowerradius = 0.0;
//       }

//       // upper radius
//       double upperradius;
//       if (idxelem < numelements - 1) {
//          upperradius = lowerradius + elementwidth;
//       } else {
//          upperradius = mymodel.UpperRadius(idx);
//       }

//       // vector to contain all nodes in element
//       auto tmp_radius = std::vector<double>(q.N(), 0.0);
//       tmp_radius[0] = lowerradius;
//       for (int idx2 = 1; idx2 < q.N() - 1; ++idx2) {
//          tmp_radius[idx2] =
//              StandardIntervalMap(q.X(idx2), lowerradius, upperradius);
//       }
//       tmp_radius[q.N() - 1] = upperradius;

//       // nodes within element
//       _vec_nodes.push_back(tmp_radius);

//       // layer in model of element
//       _vec_layers.push_back(idx);

//       if (!NotFirstElement) {
//          NotFirstElement = true;
//       }
//    }
// };

// // adding in the outer elements between the edge of the planet and the
// {
//    // get size and number of elements
//    double laydepth = ballrad - mymodel.OuterRadius();
//    assert(laydepth > 0.0 &&
//           "Ball radius must be larger than outer radius of planet");
//    int numelements = std::ceil(laydepth / max_element_size);
//    double elementwidth = laydepth / static_cast<double>(numelements);

//    // finish filling out vector
//    for (int idxelem = 0; idxelem < numelements; ++idxelem) {
//       // lower radius
//       double lowerradius = _vec_nodes.back().back();

//       // upper radius
//       double upperradius;
//       if (idxelem < numelements - 1) {
//          upperradius = lowerradius + elementwidth;
//       } else {
//          upperradius = ballrad;
//       }

//       // vector to contain all nodes in element
//       auto tmp_radius = std::vector<double>(q.N(), 0.0);
//       tmp_radius[0] = lowerradius;
//       for (int idx2 = 1; idx2 < q.N() - 1; ++idx2) {
//          tmp_radius[idx2] =
//              StandardIntervalMap(q.X(idx2), lowerradius, upperradius);
//       }
//       tmp_radius[q.N() - 1] = upperradius;

//       // nodes within element
//       _vec_nodes.push_back(tmp_radius);
//       _vec_layers.push_back(mymodel.NumberOfLayers());
//    }
// }

// all nodes
// _vec_allnodes = All_Node(_vec_noderadii, q);
// _vec_nodes = All_Node_New(_vec_noderadii, q);
// check vvec
// std::cout << "length: " << _vvec_allnodes.size()
//           << ", of individual: " << _vvec_allnodes[0].size() << "\n";
// std::cout << "nelem: " << _vec_noderadii.size() - 1 << ", npoly: " <<
// q.N()
//           << "\n";

#endif
