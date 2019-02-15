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
  unsigned int region = _eos_lg->inRegionPH(pressure, enthalpy);
  temperature = _eos_lg->temperature_from_ph(pressure, enthalpy);
  switch (region)
  {
    case 1:
      vmfrac = 0.0;
      break;

    case 2 ... 3:
      vmfrac = 1.0;
      break;

    case 4:
      {
        Real h_sat_g,h_sat_l;
        h_lat(pressure, temperature, h_sat_l, h_sat_g);
        vmfrac  = enthalpy - h_sat_l;
        vmfrac /= h_sat_g - h_sat_l;
        break;
      }

    case 5:
      vmfrac = 1.0;
      break;

    default:
      mooseError(name(), ": inRegionPH() has given an incorrect region");
  }
}

void
MoskitoPureWater2P::h_lat(
  const Real & pressure, const Real & temperature, Real & hsatl, Real & hsatg) const
{
  hsatl = _eos_lg->h_from_p_T(pressure, temperature, 1);
  hsatg = _eos_lg->h_from_p_T(pressure, temperature, 2);
}

Real
MoskitoPureWater2P::rho_g_from_p_T(Real pressure, Real temperature) const
{
  return _eos_lg->rho_from_p_T(pressure, temperature);
}

void
MoskitoPureWater2P::rho_g_from_p_T(
  Real pressure, Real temperature, Real & rho, Real & drho_dp, Real & drho_dT) const
{
  _eos_lg->rho_from_p_T(pressure, temperature, rho, drho_dp, drho_dT);
}

Real
MoskitoPureWater2P::rho_l_from_p_T(Real pressure, Real temperature) const
{
  return _eos_lg->rho_from_p_T(pressure, temperature);
}

void
MoskitoPureWater2P::rho_l_from_p_T(
  Real pressure, Real temperature, Real & rho, Real & drho_dp, Real & drho_dT) const
{
  _eos_lg->rho_from_p_T(pressure, temperature, rho, drho_dp, drho_dT);
}

Real
MoskitoPureWater2P::cp_m_from_p_T(
      const Real & pressure, const Real & temperature, const Real & vmfrac) const
{
  Real cp = 0.0;

  if (vmfrac > 0.0 && vmfrac<1.0)
  {
    cp  = _eos_lg->cp_from_p_T(pressure, temperature, 2) * vmfrac;
    cp += _eos_lg->cp_from_p_T(pressure, temperature, 1) * (1.0 - vmfrac);
  }
  else
    cp = _eos_lg->cp_from_p_T(pressure, temperature);

  return cp;
}
