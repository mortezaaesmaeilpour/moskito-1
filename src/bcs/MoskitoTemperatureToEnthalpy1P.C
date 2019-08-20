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

#include "MoskitoTemperatureToEnthalpy1P.h"

registerMooseObject("MoskitoApp", MoskitoTemperatureToEnthalpy1P);

template <>
InputParameters
validParams<MoskitoTemperatureToEnthalpy1P>()
{
  InputParameters params = validParams<NodalBC>();
  params.addRequiredParam<UserObjectName>("eos_uo",
        "The name of the userobject for EOS");
  params.addRequiredParam<Real>("temperature", "Temperature value of the BC");
  params.addRequiredCoupledVar("pressure", "Pressure nonlinear variable (Pa)");
  params.declareControllable("temperature");
  params.addClassDescription("Implements a NodalBC (Dirichlet) which calculates "
                            "specific enthalpy using temperature based on EOS "
                            " for 1 phase flow");
  return params;
}

MoskitoTemperatureToEnthalpy1P::MoskitoTemperatureToEnthalpy1P(const InputParameters & parameters)
  : NodalBC(parameters),
    _T(getParam<Real>("temperature")),
    _p(coupledValue("pressure")),
    _eos_uo(getUserObject<MoskitoEOS1P>("eos_uo"))
{
}

Real
MoskitoTemperatureToEnthalpy1P::computeQpResidual()
{
  return _u[_qp] - _eos_uo.T_to_h(_T, _p[_qp]);
}
