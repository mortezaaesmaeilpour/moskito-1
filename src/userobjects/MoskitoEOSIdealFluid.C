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
  InputParameters params = validParams<MoskitoEOS>();
  params.addParam<Real>(
      "thermal_expansion", 2.14E-4, "Constant coefficient of thermal expansion (1/K)");
  params.addParam<Real>("density0", 1000.0, "Density at zero pressure and zero temperature (kg/m^3)");
  params.addRangeCheckedParam<Real>(
      "bulk_modulus", 2.0E9, "bulk_modulus>0", "Constant bulk modulus (Pa)");

  return params;
}

MoskitoEOSIdealFluid::MoskitoEOSIdealFluid(const InputParameters & parameters)
  : MoskitoEOS(parameters),
    _thermal_expansion(getParam<Real>("thermal_expansion")),
    _bulk_modulus(getParam<Real>("bulk_modulus")),
    _density0(getParam<Real>("density0"))
{
}

Real
MoskitoEOSIdealFluid::p(Real density, Real temperature) const
{
  return _bulk_modulus * ( std::log(density / _density0)  + _thermal_expansion * temperature);
}

void
MoskitoEOSIdealFluid::dp_drhoT(
    Real density, Real temperature, Real & pressure, Real & dp_drho, Real & dp_dT) const
{
  pressure = this->p(density, temperature);
  dp_drho = _bulk_modulus / density;
  dp_dT = _bulk_modulus * _thermal_expansion;
}

void
MoskitoEOSIdealFluid::dp_drhoT_2(
    Real density, Real temperature, Real & dp_drho_2, Real & dp_dT_2) const
{
  dp_drho_2 = -_bulk_modulus / (density * density);
  dp_dT_2 = 0.0;
}
