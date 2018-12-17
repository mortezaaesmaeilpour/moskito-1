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
  : FluidProperties(parameters),
    _liquid_name(isParamValid("eos_uo_liquid") ? getParam<UserObjectName>("eos_uo_liquid")
                                           : UserObjectName(name() + ":liquid")),
    _gas_name(isParamValid("eos_uo_gas") ? getParam<UserObjectName>("eos_uo_gas")
                                         : UserObjectName(name() + ":gas"))
{
  if (!isParamValid("eos_uo_liquid"))
    if (_tid == 0 && _fe_problem.hasUserObject(_liquid_name))
      paramError("eos_uo_liquid",
                 "The two-phase fluid properties object '" + name() + "' is ",
                 "trying to create a single-phase fluid properties object with ",
                 "name '", _liquid_name, "', but a single-phase fluid properties ",
                 "object with this name already exists.");

  if (!isParamValid("eos_uo_gas"))
    if (_tid == 0 && _fe_problem.hasUserObject(_gas_name))
      paramError("eos_uo_gas",
                 "The two-phase fluid properties object '" + name() + "' is ",
                 "trying to create a single-phase fluid properties object with ",
                 "name '", _gas_name, "', but a single-phase fluid properties ",
                 "object with this name already exists.");
}
