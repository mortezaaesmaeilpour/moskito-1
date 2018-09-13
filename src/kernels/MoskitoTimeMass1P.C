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

#include "MoskitoTimeMass1P.h"

registerMooseObject("MoskitoApp", MoskitoTimeMass1P);

template <>
InputParameters
validParams<MoskitoTimeMass1P>()
{
  InputParameters params = validParams<Kernel>();

  params.addRequiredCoupledVar("enthalpy", "Specific enthalpy nonlinear variable");
  params.addClassDescription("Time derivative part of mass conservation equation for "
                  "1 phase (either liquid or gas) pipe flow and it returns pressure");

  return params;
}

MoskitoTimeMass1P::MoskitoTimeMass1P(const InputParameters & parameters)
  : Kernel(parameters),
    _h_dot(coupledDot("enthalpy")),
    _dh_dot(coupledDotDu("enthalpy")),
    _h_var_number(coupled("enthalpy")),
    _cp(getMaterialProperty<Real>("specific_heat")),
    _drho_dp(getMaterialProperty<Real>("drho_dp")),
    _drho_dp_2(getMaterialProperty<Real>("drho_dp_2")),
    _drho_dT(getMaterialProperty<Real>("drho_dT")),
    _drho_dT_2(getMaterialProperty<Real>("drho_dT_2")),
    _drho_dTdp(getMaterialProperty<Real>("drho_dTdp"))
{
}

Real
MoskitoTimeMass1P::computeQpResidual()
{
  Real r = 0.0;

  r += _drho_dp[_qp] * _u_dot[_qp];
  r += _drho_dT[_qp] * _h_dot[_qp] / _cp[_qp];
  r *= _test[_i][_qp];

  return r;
}

Real
MoskitoTimeMass1P::computeQpJacobian()
{
  Real j = 0.0;

  j += _drho_dp[_qp] * _phi[_j][_qp] * _du_dot_du[_qp];
  j += _drho_dp_2[_qp] * _phi[_j][_qp] * _u_dot[_qp];
  j += _drho_dTdp[_qp] * _phi[_j][_qp] * _h_dot[_qp] / _cp[_qp];
  j *= _test[_i][_qp];

  return j;
}

Real
MoskitoTimeMass1P::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real j = 0.0;

  if (jvar == _h_var_number)
  {
    j += _drho_dTdp[_qp] * _phi[_j][_qp] * _u_dot[_qp];
    j += _drho_dT_2[_qp] * _phi[_j][_qp] * _h_dot[_qp] / _cp[_qp];
    j += _drho_dT[_qp] * _phi[_j][_qp] * _dh_dot[_qp];
    j *= _test[_i][_qp] / _cp[_qp];
  }

  return j;
}
