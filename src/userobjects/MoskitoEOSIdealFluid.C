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

#include "MoskitoEOSIdealFluid.h"

registerMooseObject("MoskitoApp", MoskitoEOSIdealFluid);

template <>
InputParameters
validParams<MoskitoEOSIdealFluid>()
{
  InputParameters params = validParams<MoskitoEOS1P>();
  params.addParam<Real>("thermal_expansion", 4.0E-4,
        "Constant coefficient of thermal expansion (1/K)");
  params.addParam<Real>("reference_density", 998.29,
        "Density at the reference pressure and temperature (kg/m^3)");
  params.addParam<Real>("reference_temperature", 293.15,
        "Reference temperature (K)");
  params.addParam<Real>("reference_pressure", 101325,
        "Reference pressure (Pa)");
  params.addParam<Real>("reference_enthalpy", 83950,
        "Specific enthalpy (J/kg)");
  params.addParam<Real>("specific_heat", 4200,
        "Constant specific heat at constant pressure (J/kg.K)");
  params.addParam<Real>("thermal_conductivity", 0.6,
        "Constant thermal conductivity (W/m/K)");
  params.addRangeCheckedParam<Real>("bulk_modulus", 2.15E9, "bulk_modulus>0",
        "Constant bulk modulus (Pa)");

  return params;
}

MoskitoEOSIdealFluid::MoskitoEOSIdealFluid(const InputParameters & parameters)
  : MoskitoEOS1P(parameters),
    _rho_ref(getParam<Real>("reference_density")),
    _T_ref(getParam<Real>("reference_temperature")),
    _P_ref(getParam<Real>("reference_pressure")),
    _h_ref(getParam<Real>("reference_enthalpy")),
    _cp(getParam<Real>("specific_heat")),
    _lambda(getParam<Real>("thermal_conductivity")),
    _thermal_expansion(getParam<Real>("thermal_expansion")),
    _bulk_modulus(getParam<Real>("bulk_modulus"))
{
}

Real
MoskitoEOSIdealFluid::rho(const Real & pressure, const Real & temperature) const
{
  return _rho_ref * std::exp((pressure-_P_ref) / _bulk_modulus - _thermal_expansion * (temperature-_T_ref));
}

void
MoskitoEOSIdealFluid::drho_dpT(const Real & pressure, const Real & temperature,
                              Real & rho, Real & drho_dp, Real & drho_dT) const
{
  rho = this->rho(pressure, temperature);
  drho_dp = rho / _bulk_modulus;
  drho_dT = -_thermal_expansion * rho;
}

Real
MoskitoEOSIdealFluid::h_to_T(const Real & enthalpy) const
{
  return (enthalpy - _h_ref) / _cp + _T_ref;
}

Real
MoskitoEOSIdealFluid::T_to_h(const Real & temperature) const
{
  return cp(temperature) * (temperature - _T_ref) + _h_ref;
}

Real
MoskitoEOSIdealFluid::cp(const Real & temperature) const
{
  return _cp;
}

Real
MoskitoEOSIdealFluid::lambda(const Real & pressure, const Real & temperature) const
{
  return _lambda;
}
