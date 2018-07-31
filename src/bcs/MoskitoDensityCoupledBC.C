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

#include "MoskitoDensityCoupledBC.h"

registerMooseObject("MoskitoApp", MoskitoDensityCoupledBC);

template <>
InputParameters
validParams<MoskitoDensityCoupledBC>()
{
  InputParameters params = validParams<NodalBC>();
  params.addRequiredParam<FunctionName>("pressure", "Pressure nonlinear function");
  params.addRequiredParam<UserObjectName>("eos_UO", "The name of the userobject for EOS");
  params.addRequiredCoupledVar("temperature", "Temperature nonlinear variable");
  params.addClassDescription("Implements a NodalBC which calculate density based on EOS using the "
                             "coupled temperature variable and given pressure function");
  return params;
}

MoskitoDensityCoupledBC::MoskitoDensityCoupledBC(const InputParameters & parameters)
  : NodalBC(parameters),
    _T(coupledValue("temperature")),
    _T_var_number(coupled("temperature")),
    _p_func(getFunction("pressure")),
    _eos_UO(getUserObject<MoskitoEOS>("eos_UO"))
{
}

Real
MoskitoDensityCoupledBC::computeQpResidual()
{
  return _u[_qp] - _eos_UO.rho(_p_func.value(_t, *_current_node), _T[_qp]);
}

Real
MoskitoDensityCoupledBC::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real rho, drho_dp, drho_dT;
  _eos_UO.drho_dpT(_p_func.value(_t, *_current_node), _T[_qp], rho, drho_dp, drho_dT);
  if (jvar == _T_var_number)
    return -drho_dT;
  else
    return 0.0;
}
