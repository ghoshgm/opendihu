#include "spatial_discretization/finite_element_method/integrand/integrand_stiffness_matrix_laplace.h"

#include <array>

#include "utility/math_utility.h"

namespace SpatialDiscretization
{

//integrand for stiffness matrix of laplace operator, 1D
template<typename EvaluationsType,typename FunctionSpaceType,typename Term>
EvaluationsType IntegrandStiffnessMatrix<1,EvaluationsType,FunctionSpaceType,Term,Equation::hasLaplaceOperator<Term>>::
evaluateIntegrand(const Data::FiniteElements<FunctionSpaceType,Term> &data, const std::array<Vec3,1> &jacobian,
                  element_no_t elementNoLocal, const std::array<double,1> xi)
{
  EvaluationsType evaluations;

  double s = MathUtility::norm<3>(jacobian[0]);
  double integralFactor = 1. / s;

  // initialize gradient vectors of ansatz function phi_i, for node i of current element
  std::array<std::array<double,1>,FunctionSpaceType::nDofsPerElement()> gradPhi = data.functionSpace()->getGradPhi(xi);

  // loop over pairs of basis functions and evaluation integrand at xi
  for (int i=0; i<FunctionSpaceType::nDofsPerElement(); i++)
  {
    for (int j=0; j<FunctionSpaceType::nDofsPerElement(); j++)
    {
      double integrand = gradPhi[i][0] * gradPhi[j][0] * integralFactor;
      evaluations(i,j) = integrand;
      VLOG(2) << "   dofs(" <<i<< "," <<j<< ") integrand=" <<gradPhi[i][0]<< "*" <<gradPhi[j][0]<< "*" <<integralFactor<< "=" <<integrand;
    }
  }

  return evaluations;
}

//integrand for stiffness matrix of laplace operator, 2D
template<typename EvaluationsType,typename FunctionSpaceType,typename Term>
EvaluationsType IntegrandStiffnessMatrix<2,EvaluationsType,FunctionSpaceType,Term,Equation::hasLaplaceOperator<Term>>::
evaluateIntegrand(const Data::FiniteElements<FunctionSpaceType,Term> &data, const std::array<Vec3,2> &jacobian,
                  element_no_t elementNoLocal, const std::array<double,2> xi)
{
  VLOG(3) << "evaluateIntegrand laplace operator 2D";

  EvaluationsType evaluations;

  const Vec3 &zeta1 = jacobian[0];  // first column of jacobian
  const Vec3 &zetah = jacobian[1];  // second column of jacobian

  double integrationFactor = MathUtility::computeIntegrationFactor<2>(jacobian);

  double l1 = MathUtility::length<3>(zeta1);
  double lh = MathUtility::length<3>(zetah);
  double beta = acos((zeta1[0]*zetah[0] + zeta1[1]*zetah[1] + zeta1[2]*zetah[2]) / (l1 * lh));
  double alpha = MathUtility::sqr(MathUtility::PI) / (4.*beta);

  Vec3 zeta2 = cos(alpha) * zeta1 + sin(alpha) * zetah;
  double l2 = MathUtility::length<3>(zeta2);
  double l1squared = MathUtility::sqr(l1);
  double l2squared = MathUtility::sqr(l2);

  // compute the 2x2 transformation matrix T
  std::array<double,4> transformationMatrix = {
    MathUtility::sqr(cos(alpha))/l2squared + 1./l1squared,
    sin(alpha)*cos(alpha)/l2squared,
    sin(alpha)*cos(alpha)/l2squared,
    MathUtility::sqr(sin(alpha))/l2squared
  };

#ifdef DEBUG
  VLOG(3) << "transformationMatrix:";
    std::stringstream s;
  for (int i=0; i<2; i++)
  {
    for (int j=0; j<2; j++)
    {
      s << transformationMatrix[i*2+j] << " ";
    }
    s << std::endl;
  }
  VLOG(3) << s.str();
#endif

  // initialize gradient vectors of ansatz function phi_i, for node i of current element
  std::array<Vec2,FunctionSpaceType::nDofsPerElement()> gradPhi = data.functionSpace()->getGradPhi(xi);

  // loop over pairs of basis functions and evaluation integrand at xi
  for (int i=0; i<FunctionSpaceType::nDofsPerElement(); i++)
  {
    for (int j=0; j<FunctionSpaceType::nDofsPerElement(); j++)
    {
      double integrand = MathUtility::applyTransformation(transformationMatrix, gradPhi[i], gradPhi[j]) * integrationFactor;
      evaluations(i,j) = integrand;
    }
  }

  return evaluations;
}

//integrand for stiffness matrix of laplace operator, 3D
template<typename EvaluationsType,typename FunctionSpaceType,typename Term>
EvaluationsType IntegrandStiffnessMatrix<3,EvaluationsType,FunctionSpaceType,Term,Equation::hasLaplaceOperator<Term>>::
evaluateIntegrand(const Data::FiniteElements<FunctionSpaceType,Term> &data, const std::array<Vec3,3> &jacobian,
                  element_no_t elementNoLocal, const std::array<double,3> xi)
{
  EvaluationsType evaluations;

  // compute the 3x3 transformation matrix T = J^{-1}J^{-T} and the absolute of the determinant of the jacobian
  double determinant;
  std::array<double,9> transformationMatrix = MathUtility::computeTransformationMatrixAndDeterminant(jacobian, determinant);

#ifdef DEBUG
  VLOG(3) << "transformationMatrix:";
  std::stringstream s;
  for (int i=0; i<3; i++)
  {
    for (int j=0; j<3; j++)
    {
      s << transformationMatrix[i*3+j] << " ";
    }
    s << std::endl;
  }
  VLOG(3) << s.str();
#endif

  const std::array<Vec3,FunctionSpaceType::nDofsPerElement()> gradPhi = data.functionSpace()->getGradPhi(xi);

  // loop over pairs of basis functions and evaluation integrand at xi
  for (int i=0; i<FunctionSpaceType::nDofsPerElement(); i++)
  {
    for (int j=0; j<FunctionSpaceType::nDofsPerElement(); j++)
    {

      //! computes gradPhi[i]^T * T * gradPhi[j] where T is the symmetric transformation matrix
      double integrand = MathUtility::applyTransformation(transformationMatrix, gradPhi[i], gradPhi[j]) * fabs(determinant);
      evaluations(i,j) = integrand;
    }
  }

  return evaluations;
}

} // namespace
