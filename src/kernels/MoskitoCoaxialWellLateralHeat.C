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

#include "MoskitoCoaxialWellLateralHeat.h"

registerMooseObject("MooseApp", MoskitoCoaxialWellLateralHeat);

template <>
InputParameters
validParams<MoskitoCoaxialWellLateralHeat>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Source term based on the Fourier Law for coaxial "
                             "well for injection/production in the inner and "
                             "outer string.");
  return params;
}

MoskitoCoaxialWellLateralHeat::MoskitoCoaxialWellLateralHeat(const InputParameters & parameters)
  : Kernel(parameters),
{
}

Real
MoskitoCoaxialWellLateralHeat::computeQpResidual()
{
  Real factor = _scale * _postprocessor * _function.value(_t, _q_point[_qp]);
  return _test[_i][_qp] * -factor;
}
