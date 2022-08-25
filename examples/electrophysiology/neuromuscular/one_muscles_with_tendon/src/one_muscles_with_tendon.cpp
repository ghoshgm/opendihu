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

// Fast monodomain: 0D / 1D
using MonodomainSolver =
  FastMonodomainSolver<                               // a wrapper that improves performance of multidomain
    Control::MultipleInstances<                       // fibers
      OperatorSplitting::Strang<
        Control::MultipleInstances<
          TimeSteppingScheme::Heun<                   // fiber reaction term
            CellmlAdapter<
              N_STATES, N_ALGEBRAICS,                 // depends on the cellml model
              FunctionSpace::FunctionSpace<
                Mesh::StructuredDeformableOfDimension<1>,
                BasisFunction::LagrangeOfOrder<1>
              >
            >
          >
        >,
        Control::MultipleInstances<
          FiberDiffusionSolver<                       // fiber diffusion, note that implicit euler gives lower error in this case than crank nicolson
            SpatialDiscretization::FiniteElementMethod<
              Mesh::StructuredDeformableOfDimension<1>,
              BasisFunction::LagrangeOfOrder<1>,
              Quadrature::Gauss<2>,
              Equation::Dynamic::IsotropicDiffusion
            >
          >
        >
      >
    >
  >;

// 3D bidomain
using BidomainSolver =
  OutputWriter::OutputSurface<
    TimeSteppingScheme::StaticBidomainSolver<         // bidomain
      SpatialDiscretization::FiniteElementMethod<     // FEM for initial potential flow, fiber directions
        Mesh::StructuredDeformableOfDimension<3>,
        BasisFunction::LagrangeOfOrder<1>,
        Quadrature::Gauss<3>,
        Equation::Static::Laplace
      >,
      SpatialDiscretization::FiniteElementMethod<     // anisotropic diffusion
        Mesh::StructuredDeformableOfDimension<3>,
        BasisFunction::LagrangeOfOrder<1>,
        Quadrature::Gauss<5>,
        Equation::Dynamic::DirectionalDiffusion
      >
    >
  >;

// define helper function space for various activation signals, this is actually a vector space
typedef FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<1>,BasisFunction::LagrangeOfOrder<1>> HelperFunctionSpace;



int main(int argc, char *argv[])
{
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);

  // Hierarchical coupling to use small time steps for spindles,
  // motoneurons, monodomain, and bidomain coupled to the contraction solver
  // with large time steps
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
      // motoneuron solver
      TimeSteppingScheme::Heun<
        CellmlAdapter<
          6,14,  // nStates,nAlgebraics
          HelperFunctionSpace
        >
      >,
      // map from λ in the 3D mesh to muscle spindles input
      Control::MapDofs<
        HelperFunctionSpace,
      
        // map from motoneuronMesh to stimulated nodes
        Control::MapDofs<
          HelperFunctionSpace,
                  
          // electrophysiology + mechanics
          Control::Coupling<

            // Fast monodomain: 0D / 1D + 3D bidomain
            Control::MapDofs< // motor neuron -> activation signal
                HelperFunctionSpace,
                Control::Coupling<
                    MonodomainSolver,
                    BidomainSolver
                >
            >,   
            // muscle contraction
            MuscleContractionSolver<
                Mesh::StructuredDeformableOfDimension<3>
            >    
        >
       >
      >
     >
  problem(settings);
  
  problem.run();
  
  return EXIT_SUCCESS;
}
