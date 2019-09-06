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

#include "MoskitoPureWater1P.h"

registerMooseObject("MoskitoApp", MoskitoPureWater1P);

template <>
InputParameters
validParams<MoskitoPureWater1P>()
{
  InputParameters params = validParams<MoskitoEOS1P>();
  params.addParam<Real>("thermal_conductivity", 0.6,
        "Constant thermal conductivity (W/m/K)");

  return params;
}

MoskitoPureWater1P::MoskitoPureWater1P(const InputParameters & parameters)
  : MoskitoEOS1P(parameters),
    _lambda(getParam<Real>("thermal_conductivity"))
{
  std::string eos_name = UserObjectName(name() + ":LiquidGas");
  {
    std::string class_name = "MoskitoWater97FluidProperties";
    InputParameters params = _app.getFactory().getValidParams(class_name);
    _fe_problem.addUserObject(class_name, eos_name, params);
  }
  _eos_1P = &_fe_problem.getUserObjectTempl<MoskitoWater97FluidProperties>(eos_name);
}

Real
MoskitoPureWater1P::rho_from_p_T(const Real & pressure, const Real & temperature) const
{
  Real rho = _eos_1P->rho_from_p_T(pressure, temperature);
  return rho;
}

void
MoskitoPureWater1P::rho_from_p_T(const Real & pressure, const Real & temperature,
                              Real & rho, Real & drho_dp, Real & drho_dT) const
{
      _eos_1P->rho_from_p_T(pressure, temperature, rho, drho_dp, drho_dT);
}

Real
MoskitoPureWater1P::h_to_T(const Real & pressure, const Real & enthalpy) const
{
  unsigned int region = _eos_1P->inRegionPH(pressure, enthalpy);

  Real temp;

  if (region == 4)
    mooseError(name(), "This EOS is not valid for 2 phase region");

  else
    temp = _eos_1P->temperature_from_ph(pressure, enthalpy);

  return temp;
}

Real
MoskitoPureWater1P::T_to_h(const Real & pressure, const Real & temperature) const
{
  Real enthalpy = _eos_1P->h_from_p_T(pressure, temperature);

  return enthalpy;
}

Real
MoskitoPureWater1P::cp(const Real & pressure, const Real & temperature) const
{
  Real cp;
  cp = _eos_1P->cp_from_p_T(pressure, temperature);
  return _eos_1P->cp_from_p_T(pressure, temperature);
}

Real
MoskitoPureWater1P::lambda(const Real & pressure, const Real & temperature) const
{
return _lambda;
}
