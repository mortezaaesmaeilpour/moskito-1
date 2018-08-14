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
  params.addParam<Real>("molar_mass", 29.0e-3, "Constant molar mass of the fluid (kg/mol)");
  params.addRequiredParam<Real>("R", "Gas constant");

  return params;
}

MoskitoEOSIdealGas::MoskitoEOSIdealGas(const InputParameters & parameters)
  : MoskitoEOS(parameters), _molar_mass(getParam<Real>("molar_mass")), _R(getParam<Real>("R"))
{
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