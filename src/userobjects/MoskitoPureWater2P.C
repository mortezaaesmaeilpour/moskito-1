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

#include "MoskitoPureWater2P.h"

registerMooseObject("MoskitoApp", MoskitoPureWater2P);

template <>
InputParameters
validParams<MoskitoPureWater2P>()
{
  InputParameters params = validParams<MoskitoEOS2P>();
  return params;
}

MoskitoPureWater2P::MoskitoPureWater2P(const InputParameters & parameters)
  : MoskitoEOS2P(parameters)
{
  std::string eos_name = UserObjectName(name() + ":LiquidGas");
  {
    std::string class_name = "MoskitoWater97FluidProperties";
    InputParameters params = _app.getFactory().getValidParams(class_name);
    _fe_problem.addUserObject(class_name, eos_name, params);
  }
  _eos_lg = &_fe_problem.getUserObject<MoskitoWater97FluidProperties>(eos_name);
}

void
MoskitoPureWater2P::VMFrac_from_p_h(
  const Real & pressure, const Real & enthalpy, Real & vmfrac, Real & temperature) const
{

}

void
MoskitoPureWater2P::h_lat(
  const Real & pressure, const Real & temperature, Real & hsatl, Real & hsatg) const
{
  
}

Real
MoskitoPureWater2P::rho_g_from_p_T(Real pressure, Real temperature) const
{
  return 0;
}

void
MoskitoPureWater2P::rho_g_from_p_T(
  Real pressure, Real temperature, Real & rho, Real & drho_dp, Real & drho_dT) const
{

}

Real
MoskitoPureWater2P::rho_l_from_p_T(Real pressure, Real temperature) const
{
  return 0;
}

void
MoskitoPureWater2P::rho_l_from_p_T(
  Real pressure, Real temperature, Real & rho, Real & drho_dp, Real & drho_dT) const
{

}

Real
MoskitoPureWater2P::cp_m_from_p_T(
      const Real & pressure, const Real & temperature, const Real & vmfrac) const
{
  return 0;
}
