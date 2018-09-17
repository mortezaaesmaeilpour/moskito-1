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
  InputParameters params = validParams<MoskitoEOS>();

  params.addRequiredParam<Real>("molar_mass", "Molar mass of the gas (kg/mol)");
  params.addParam<Real>(
      "cp", 1.005e3, "Constant specific heat capacity at constant pressure (J/kg/K)");

  return params;
}

MoskitoEOSIdealGas::MoskitoEOSIdealGas(const InputParameters & parameters)
  : MoskitoEOS(parameters), _molar_mass(getParam<Real>("molar_mass")), _R(8.3144598)
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

void
MoskitoEOSIdealGas::drho_dpT_2(
    Real pressure, Real temperature, Real & drho_dp_2, Real & drho_dT_2, Real & drho_dTdp) const
{
  Real rho = this->rho(pressure, temperature);
  drho_dp_2 = 0.0;
  drho_dT_2 = 2.0 * pressure * _molar_mass / (_R * temperature * temperature * temperature);
  drho_dTdp = -_molar_mass / (_R * temperature * temperature);
}

Real
MoskitoEOSIdealGas::p(Real density, Real temperature) const
{
  return density * _R * temperature / _molar_mass;
}

void
MoskitoEOSIdealGas::dp_drhoT(
    Real density, Real temperature, Real & pressure, Real & dp_drho, Real & dp_dT) const
{
  pressure = this->p(density, temperature);
  dp_drho = _R * temperature / _molar_mass;
  dp_dT = density * _R / _molar_mass;
}

void
MoskitoEOSIdealGas::dp_drhoT_2(Real density,
                               Real temperature,
                               Real & dp_drho_2,
                               Real & dp_dT_2) const
{
  dp_drho_2 = 0.0;
  dp_dT_2 = 0.0;
}

Real
MoskitoEOSIdealGas::h_to_T(Real enthalpy) const
{
  return enthalpy / _cp;
}
