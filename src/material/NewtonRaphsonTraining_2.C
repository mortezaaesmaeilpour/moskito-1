#include "NewtonRaphsonTraining_2.h"

registerMooseObject("MoskitoApp", NewtonRaphsonTraining_2);

template <>
InputParameters
validParams<NewtonRaphsonTraining_2>()
{
    InputParameters params = validParams<Material>();
    params += validParams<NewtonIteration>();
    return params;
}

NewtonRaphsonTraining_2::NewtonRaphsonTraining_2(const InputParameters & parameters)
  : Material(parameters),
    NewtonIteration(parameters)
  {

  }
Real
NewtonRaphsonTraining_2::calculate_Reynolds(Real rho, Real v, Real kin_vis, Real length)
{
  Real Reynolds;
  Reynolds = rho * v * length / kin_vis;
  return Reynolds;
}

Real
NewtonRaphsonTraining_2::initialGuess(const Real trial_value)
{
  return 100.0;
}

// Real
// NewtonIteration::maximumPermissibleValue(const Real trial_value) const
// {
//   return 1000.0;
// }

Real
NewtonRaphsonTraining_2::minimumPermissibleValue(const Real trail_value) const
{
    return 0.00000001;
}

Real
NewtonRaphsonTraining_2::computeReferenceResidual(const Real trial_value, const Real scalar)
{
  return 0.001;
}

Real
NewtonRaphsonTraining_2::computeResidual(const Real trial_value, const Real scalar)
{

  Real epsi = 0.001;
  Real dia = 0.1;

  Real rho = 1000.0; // kg/m³
  Real v = 3.0; // m/s
  Real kin_vis = 1.0; // mPa s;
  Real length = 10.0; // m

  Real Re = calculate_Reynolds(rho, v, kin_vis, length);

  Real bracket = 0.0, r = 0.0;
  bracket = epsi / (3.7 * dia);
  bracket += 2.51 / (Re * std::pow(scalar,0.5));

  r =  std::log10(bracket) * ( - 2.0);
  r -= 1 / std::pow(scalar,0.5);
  return r;
}

Real
NewtonRaphsonTraining_2::computeDerivative(const Real trial_value, const Real scalar)
{
  Real bracket = 0.0, f_dash = 0.0;

  Real epsi = 0.001;
  Real dia = 0.1;

  Real rho = 1000.0; // kg/m³
  Real v = 3.0; // m/s
  Real kin_vis = 1.0; // mPa s;
  Real length = 10.0; // m

  Real Re = calculate_Reynolds(rho, v, kin_vis, length);

  bracket = epsi / (3.7 * dia);
  bracket += 2.51 / (Re * std::pow(scalar,0.5));

  f_dash = ( - 0.5) * 2.5 / (Re * std::pow(scalar, 3.0 / 2.0));
  f_dash *= - 2.0 / (bracket * std::log(10));
  f_dash += 1.0 / (2.0 * std::pow(scalar, 3.0 / 2.0));
  return f_dash;
}

void
NewtonRaphsonTraining_2::computeQpProperties()
{
  Real rho = 1000.0; // kg/m³
  Real v = 3.0; // m/s
  Real kin_vis = 1.0; // mPa s;
  Real length = 10.0; // m

  Real Re = calculate_Reynolds(rho, v, kin_vis, length);

  Real epsi = 0.1;
  Real dia = 10;
  Real Rel_rough = epsi / dia;

  Real scalar, trial = 0.0;
  returnNewtonSolve(trial, scalar, _console);

}
