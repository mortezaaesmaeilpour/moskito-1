/**************************************************************************/
/*  MOSKITO - Multiphysics cOupled Simulator toolKIT for wellbOres        */
/*                                                                        */
/*  Copyright (C) 2017 by Maziar Gholami Korzani                          */
/*  Karlsruhe Institute of Technology, Institute of Applied Geosciences   */
/*  Division of Geothermal Research                                       */
/*                                                                        */
/*  This file is part of MOSKITO App                                      */
/*                                                                        */
/*  This program is free software: you can redistribute it and/or modify  */
/*  it under the terms of the GNU General Public License as published by  */
/*  the Free Software Foundation, either version 3 of the License, or     */
/*  (at your option) any later version.                                   */
/*                                                                        */
/*  This program is distributed in the hope that it will be useful,       */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          */
/*  GNU General Public License for more details.                          */
/*                                                                        */
/*  You should have received a copy of the GNU General Public License     */
/*  along with this program.  If not, see <http://www.gnu.org/licenses/>  */
/**************************************************************************/

#include "MoskitoTemperatureToEnthalpy2P.h"

registerMooseObject("MoskitoApp", MoskitoTemperatureToEnthalpy2P);

template <>
InputParameters
validParams<MoskitoTemperatureToEnthalpy2P>()
{
  InputParameters params = validParams<NodalBC>();
  // params.addRequiredParam<UserObjectName>("eos_uo",
  //       "The name of the userobject for EOS");
  params.addRequiredParam<Real>("temperature", "Temperature value of the BC");
  params.addRequiredParam<Real>("mass_fraction", "Gas to Liquid ratio; value between 0 and 1"
                                "0.0 = Liquid phase + critical; 1.0 = vapour phase ");
  params.addRequiredCoupledVar("pressure", "Pressure nonlinear variable (Pa)");
  params.declareControllable("temperature");
  params.addClassDescription("Implements a NodalBC (Dirichlet) which calculates "
                            "specific enthalpy using temperature based on EOS "
                            " for 2 phase flow in Pure water");
  return params;
}

MoskitoTemperatureToEnthalpy2P::MoskitoTemperatureToEnthalpy2P(const InputParameters & parameters)
  : NodalBC(parameters),
    _T(getParam<Real>("temperature")),
    _vmfrac(getParam<Real>("mass_fraction")),
    _p(coupledValue("pressure"))
{
  std::string eos_name = UserObjectName(name() + ":LiquidGas");
  {
    std::string class_name = "MoskitoWater97FluidProperties";
    InputParameters params = _app.getFactory().getValidParams(class_name);
    _fe_problem.addUserObject(class_name, eos_name, params);
  }
  _eos_uo = &_fe_problem.getUserObjectTempl<MoskitoWater97FluidProperties>(eos_name);
}

Real
MoskitoTemperatureToEnthalpy2P::computeQpResidual()
{
  Real enthalpy = 0.0;
  if (_vmfrac != 1.0 &&  _vmfrac != 0.0)
  {
    Real pressure = _eos_uo->vaporPressure(_T);
    Real hsatg = _eos_uo->h_from_p_T(pressure, _T, 2);
    Real hsatl = _eos_uo->h_from_p_T(pressure, _T, 1);
    enthalpy = hsatl + (hsatg - hsatl) * _vmfrac;
  }

  else
  {
    unsigned int region = _eos_uo->inRegion(_p[_qp], _T);
    enthalpy = _eos_uo->h_from_p_T(_p[_qp], _T, region);
  }
  return _u[_qp] - enthalpy;
}
