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

#include "MoskitoEOS1P_PureWater.h"

registerMooseObject("MoskitoApp", MoskitoEOS1P_PureWater);

template <>
InputParameters
validParams<MoskitoEOS1P_PureWater>()
{
  InputParameters params = validParams<MoskitoEOS1P>();

  return params;
}

MoskitoEOS1P_PureWater::MoskitoEOS1P_PureWater(const InputParameters & parameters)
  : MoskitoEOS1P(parameters)
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
MoskitoEOS1P_PureWater::rho_from_p_T(const Real & pressure, const Real & temperature) const
{
  return _eos_1P->rho_from_p_T(pressure, temperature);
}

void
MoskitoEOS1P_PureWater::rho_from_p_T(const Real & pressure, const Real & temperature,
                              Real & rho, Real & drho_dp, Real & drho_dT) const
{
  _eos_1P->rho_from_p_T(pressure, temperature, rho, drho_dp, drho_dT);
}

Real
MoskitoEOS1P_PureWater::h_to_T(const Real & pressure, const Real & enthalpy) const
{
  unsigned int region = _eos_1P->inRegionPH(pressure, enthalpy);
  if (region == 4)
    mooseError(name(), "This EOS is not valid for 2 phase region");

  return _eos_1P->temperature_from_ph(pressure, enthalpy);
}

Real
MoskitoEOS1P_PureWater::T_to_h(const Real & pressure, const Real & temperature) const
{
  return _eos_1P->h_from_p_T(pressure, temperature);
}

Real
MoskitoEOS1P_PureWater::cp(const Real & pressure, const Real & temperature) const
{
  return _eos_1P->cp_from_p_T(pressure, temperature);
}

Real
MoskitoEOS1P_PureWater::lambda(const Real & pressure, const Real & temperature) const
{
return _eos_1P->k_from_p_T(pressure, temperature);
}
