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

#include "MoskitoEOS2P.h"

registerMooseObject("MoskitoApp", MoskitoEOS2P);

template <>
InputParameters
validParams<MoskitoEOS2P>()
{
  InputParameters params = validParams<GeneralUserObject>();
  params.addRequiredParam<UserObjectName>("eos_uo_gas",
        "The name of the EOS userobject for gas");
  params.addRequiredParam<UserObjectName>("eos_uo_liquid",
        "The name of the EOS userobject for liquid");
  return params;
}

MoskitoEOS2P::MoskitoEOS2P(const InputParameters & parameters)
  : GeneralUserObject(parameters),
    gas(getUserObject<MoskitoEOS1P>("eos_uo_gas")),
    liquid(getUserObject<MoskitoEOS1P>("eos_uo_liquid"))
{
}

MoskitoEOS2P::~MoskitoEOS2P() {}

Real
MoskitoEOS2P::GasMassFraction(const Real & enthalpy, const Real & pressure) const
{
  return 0.0;
}

Real
MoskitoEOS2P::cp(const Real & massfraction, const Real & temperature) const
{
  return gas.cp(temperature) * massfraction + liquid.cp(temperature) * (1.0 - massfraction);
}

Real
MoskitoEOS2P::h_to_T(const Real & enthalpy, const Real & pressure) const
{
    return 0.0;
}
