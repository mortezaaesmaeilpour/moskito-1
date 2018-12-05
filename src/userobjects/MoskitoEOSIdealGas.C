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

#include "MoskitoEOSIdealGas.h"

registerMooseObject("MoskitoApp", MoskitoEOSIdealGas);

template <>
InputParameters
validParams<MoskitoEOSIdealGas>()
{
  InputParameters params = validParams<MoskitoEOS1P>();

  params.addRequiredParam<Real>("molar_mass", "Molar mass of the gas (kg/mol)");
  params.addParam<Real>(
      "cp", 1.005e3, "Constant specific heat capacity at constant pressure (J/kg/K)");

  return params;
}

MoskitoEOSIdealGas::MoskitoEOSIdealGas(const InputParameters & parameters)
  : MoskitoEOS1P(parameters), _molar_mass(getParam<Real>("molar_mass")), _R(8.3144598)
{
  _cp = getParam<Real>("cp");
  _density_ref = 0.0;
  _T_ref = 0.0;
  _P_ref = 0.0;
  _h_ref = 0.0;
}

Real
MoskitoEOSIdealGas::rho(Real pressure, Real temperature) const
{
  return pressure * _molar_mass / (_R * temperature);
}

void
MoskitoEOSIdealGas::drho_dpT(
    Real pressure, Real temperature, Real & rho, Real & drho_dp, Real & drho_dT) const
{
  rho = this->rho(pressure, temperature);
  drho_dp = _molar_mass / (_R * temperature);
  drho_dT = -pressure * _molar_mass / (_R * temperature * temperature);
}

Real
MoskitoEOSIdealGas::h_to_T(Real enthalpy) const
{
  return enthalpy / _cp;
}

Real
MoskitoEOSIdealGas::T_to_h(Real temperature) const
{
  return _cp * temperature;
}
