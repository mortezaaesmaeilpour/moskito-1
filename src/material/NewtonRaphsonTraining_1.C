#include "NewtonRaphsonTraining_1.h"

registerMooseObject("MoskitoApp", NewtonRaphsonTraining_1);

template <>
InputParameters
validParams<NewtonRaphsonTraining_1>()
{
    InputParameters params = validParams<Material>();
    params += validParams<NewtonIteration>();
    return params;
}

NewtonRaphsonTraining_1::NewtonRaphsonTraining_1(const InputParameters & parameters)
  : Material(parameters),
    NewtonIteration(parameters)
  {

  }

Real
NewtonRaphsonTraining_1::initialGuess(const Real trial_value)
{
  return 1.0;
}

Real
NewtonRaphsonTraining_1::computeReferenceResidual(const Real trial_value, const Real scalar)
{
  return 1.0;
}

Real
NewtonRaphsonTraining_1::computeResidual(const Real trial_value, const Real scalar)
{
  Real r = 0.0;
  r = - scalar * scalar * scalar;
  r -= 4.0 * scalar;
  r += 10.0;
  return r;
}

Real
NewtonRaphsonTraining_1::computeDerivative(const Real trial_value, const Real scalar)
{
  Real f_dash = 0.0;
  f_dash = - 3.0 *scalar *scalar;
  f_dash -= 4.0;
  return f_dash;
}

void
NewtonRaphsonTraining_1::computeQpProperties()
{
  Real scalar, trial = 0.0;
  returnNewtonSolve(trial, scalar, _console);
  // std::cout<<"z coord ="<<_q_point[_qp] <<std::endl;
  std::cout<<"x0 = "<< scalar <<std::endl;
}
