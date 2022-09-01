#include <iostream>
#include <cstdlib>

#include "opendihu.h"



////// Possible options

// #define Shorten
// #define HodgkinHuxley
#define HodgkinHuxlexRazumova

// #define FiberDiffusionSolver TimeSteppingScheme::ImplicitEuler
// #define FiberDiffusionSolver TimeSteppingScheme::CrankNicolson
template<typename T> using FiberDiffusionSolver = TimeSteppingScheme::ImplicitEuler<T>;
// template<typename T> using FiberDiffusionSolver = TimeSteppingScheme::CrankNicolson<T>;

////// ^^^^^^^^^^^^^^^^

#ifdef Shorten
#define N_STATES 57
#define N_ALGEBRAICS 1
#endif
#ifdef HodgkinHuxley
#define N_STATES 4
#define N_ALGEBRAICS 9
#endif
#ifdef HodgkinHuxlexRazumova
#define N_STATES 9
#define N_ALGEBRAICS 19
#endif

typedef Mesh::StructuredDeformableOfDimension<3> MeshType;

// define helper function space for various activation signals, this is actually a vector space
typedef FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<1>,BasisFunction::LagrangeOfOrder<1>> HelperFunctionSpace;



int main(int argc, char *argv[])
{
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);

  // Hierarchical coupling to use small time steps for spindles, golgi,
  // motoneurons, monodomain, and bidomain coupled to the contraction solver
  // with large time steps
  Control::Coupling<
    Control::MultipleCoupling<
      // muscle spindles solver
      Control::MapDofs< // spindles -> motor neurons
        HelperFunctionSpace,
        TimeSteppingScheme::Heun<
          CellmlAdapter<
            8,23,
            HelperFunctionSpace
          >
        >
      >,
      // golgi tendon organ
      Control::MapDofs< // golgi tendon organs -> interneurons
        HelperFunctionSpace,
        TimeSteppingScheme::Heun<
          CellmlAdapter<
            4,9,
            HelperFunctionSpace
          >
        >
      >,
      // interneurons
      Control::MapDofs< // interneurons -> motor neurons
        HelperFunctionSpace,
        TimeSteppingScheme::Heun<
          CellmlAdapter<
            4,9,
            HelperFunctionSpace
          >
        >
      >,
      // motor neurons
      Control::MapDofs< // signals from spindles and interneurons -> motor neuron input
        HelperFunctionSpace,
        TimeSteppingScheme::Heun<
          CellmlAdapter<
            6,14,
            HelperFunctionSpace
          >
        >
      >,
      // Fast monodomain: 0D / 1D + 3D bidomain
      Control::MapDofs< // motor neuron -> activation signal
        HelperFunctionSpace,
            // Multidomain, Strang splitting of CellmlAdapter and MultidomainWithFatSolver
            OperatorSplitting::Strang<
              Control::MultipleInstances<
                TimeSteppingScheme::Heun<
                  CellmlAdapter<
                    N_STATES, N_ALGEBRAICS,                 // depends on the cellml model
                    FunctionSpace::FunctionSpace<MeshType,BasisFunction::LagrangeOfOrder<1>>  // same function space as for anisotropic diffusion
                  >  
                >
              >,
              OutputWriter::OutputSurface<
                TimeSteppingScheme::MultidomainWithFatSolver<       // multidomain
                  SpatialDiscretization::FiniteElementMethod<       //FEM for initial potential flow, fibre directions
                    MeshType,
                    BasisFunction::LagrangeOfOrder<1>,
                    Quadrature::Gauss<3>,
                    Equation::Static::Laplace
                  >,
                  SpatialDiscretization::FiniteElementMethod<       // anisotropic diffusion
                    MeshType,
                    BasisFunction::LagrangeOfOrder<1>,
                    Quadrature::Gauss<3>,
                    Equation::Dynamic::DirectionalDiffusion
                  >,
                  SpatialDiscretization::FiniteElementMethod<       // isotropic diffusion in fat layer
                    MeshType,
                    BasisFunction::LagrangeOfOrder<1>,
                    Quadrature::Gauss<3>,
                    Equation::Dynamic::IsotropicDiffusion
                  >
                >
              >
            >
      >,
      // Fast monodomain: 0D / 1D + 3D bidomain
      Control::MapDofs< // motor neuron -> activation signal
        HelperFunctionSpace,
            // Multidomain, Strang splitting of CellmlAdapter and MultidomainWithFatSolver
            OperatorSplitting::Strang<
              Control::MultipleInstances<
                TimeSteppingScheme::Heun<
                  CellmlAdapter<
                    N_STATES, N_ALGEBRAICS,                 // depends on the cellml model
                    FunctionSpace::FunctionSpace<MeshType,BasisFunction::LagrangeOfOrder<1>>  // same function space as for anisotropic diffusion
                  >  
                >
              >,
              OutputWriter::OutputSurface<
                TimeSteppingScheme::MultidomainWithFatSolver<       // multidomain
                  SpatialDiscretization::FiniteElementMethod<       //FEM for initial potential flow, fibre directions
                    MeshType,
                    BasisFunction::LagrangeOfOrder<1>,
                    Quadrature::Gauss<3>,
                    Equation::Static::Laplace
                  >,
                  SpatialDiscretization::FiniteElementMethod<       // anisotropic diffusion
                    MeshType,
                    BasisFunction::LagrangeOfOrder<1>,
                    Quadrature::Gauss<3>,
                    Equation::Dynamic::DirectionalDiffusion
                  >,
                  SpatialDiscretization::FiniteElementMethod<       // isotropic diffusion in fat layer
                    MeshType,
                    BasisFunction::LagrangeOfOrder<1>,
                    Quadrature::Gauss<3>,
                    Equation::Dynamic::IsotropicDiffusion
                  >
                >
              >
            >
      >
    >,
    // 2x mechanics TODO make sure that no data is copied
    Control::MapDofs< // lambda -> muscle splindle input
      HelperFunctionSpace,
      Control::MapDofs< // T -> golgi tendon organ input
        HelperFunctionSpace,
        Control::Coupling<
          MuscleContractionSolver<
            Mesh::StructuredDeformableOfDimension<3>
          >,
          MuscleContractionSolver<
            Mesh::StructuredDeformableOfDimension<3>
          >
        >
      >
    >
  > problem(settings);
  
  problem.run();
  
  return EXIT_SUCCESS;
}
