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

#include "MoskitoTimeEnergy.h"

registerMooseObject("MoskitoApp", MoskitoTimeEnergy);

template <>
InputParameters
validParams<MoskitoTimeEnergy>()
{
  InputParameters params = validParams<Kernel>();

  params.addRequiredCoupledVar("pressure", "Pressure nonlinear variable");
  params.addRequiredCoupledVar("flowrate", "Volumetric flowrate nonlinear variable");
  params.addClassDescription("Time derivative part of energy conservation equation for "
                  "1 phase (either liquid or gas) pipe flow and it returns enthalpy");

  return params;
}

MoskitoTimeEnergy::MoskitoTimeEnergy(const InputParameters & parameters)
  : Kernel(parameters),
    _q(coupledValue("flowrate")),
    _p_dot(coupledDot("pressure")),
    _q_dot(coupledDot("flowrate")),
    _dp_dot(coupledDotDu("pressure")),
    _dq_dot(coupledDotDu("flowrate")),
    _p_var_number(coupled("pressure")),
    _q_var_number(coupled("flowrate")),
    _area(getMaterialProperty<Real>("well_area")),
    _cp(getMaterialProperty<Real>("specific_heat")),
    _rho(getMaterialProperty<Real>("density")),
    _drho_dp(getMaterialProperty<Real>("drho_dp")),
    _drho_dp_2(getMaterialProperty<Real>("drho_dp_2")),
    _drho_dT(getMaterialProperty<Real>("drho_dT")),
    _drho_dT_2(getMaterialProperty<Real>("drho_dT_2")),
    _drho_dTdp(getMaterialProperty<Real>("drho_dTdp"))
{
}

Real
MoskitoTimeEnergy::computeQpResidual()
{
  Real r = 0.0;

  r += _drho_dp[_qp] * _p_dot[_qp];
  r += _drho_dT[_qp] * _u_dot[_qp] / _cp[_qp];
  r *= (_u[_qp] + _q[_qp] * _q[_qp] / (2.0 * _area[_qp] * _area[_qp]));
  r += _rho[_qp] * _u_dot[_qp];
  r += _rho[_qp] * _q[_qp] * _q_dot[_qp] / (_area[_qp] * _area[_qp]);
  r -= _p_dot[_qp];
  r *= _test[_i][_qp];

  return r;
}

Real
MoskitoTimeEnergy::computeQpJacobian()
{
  Real j = 0.0;

  j += _drho_dTdp[_qp] * _phi[_j][_qp] * _p_dot[_qp] / _cp[_qp];
  j += _drho_dT_2[_qp] * _phi[_j][_qp] * _u_dot[_qp] / (_cp[_qp] * _cp[_qp]);
  j += _drho_dT[_qp] * _phi[_j][_qp] * _du_dot_du[_qp] / _cp[_qp];
  j *= (_u[_qp] + _q[_qp] * _q[_qp] / (2.0 * _area[_qp] * _area[_qp]));
  j += (_drho_dp[_qp] * _p_dot[_qp] + _drho_dT[_qp] * _u_dot[_qp] / _cp[_qp]) * _phi[_j][_qp];
  j += _drho_dT[_qp] * _phi[_j][_qp] * _u_dot[_qp] / _cp[_qp];
  j += _rho[_qp] * _phi[_j][_qp] * _du_dot_du[_qp];
  j += _drho_dT[_qp] * _phi[_j][_qp] * _q[_qp] * _q_dot[_qp] / (_area[_qp] * _area[_qp] * _cp[_qp]);
  j *= _test[_i][_qp];

  return j;
}

Real
MoskitoTimeEnergy::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real j = 0.0;

  if (jvar == _q_var_number)
  {
    j += _drho_dp[_qp] * _p_dot[_qp];
    j += _drho_dT[_qp] * _u_dot[_qp] / _cp[_qp];
    j *= _q[_qp] * _phi[_j][_qp] / (_area[_qp] * _area[_qp]);
    j += _rho[_qp] * _phi[_j][_qp] * _q_dot[_qp] / (_area[_qp] * _area[_qp]);
    j += _rho[_qp] * _q[_qp] * _phi[_j][_qp] * _dq_dot[_qp] / (_area[_qp] * _area[_qp]);
    j *= _test[_i][_qp];
  }

  if (jvar == _p_var_number)
  {
    j += _drho_dp_2[_qp] * _phi[_j][_qp] * _p_dot[_qp];
    j += _drho_dp[_qp] * _phi[_j][_qp] * _dp_dot[_qp];
    j += _drho_dTdp[_qp] * _phi[_j][_qp] * _u_dot[_qp] / _cp[_qp];
    j *= (_u[_qp] + _q[_qp] * _q[_qp] / (2.0 * _area[_qp] * _area[_qp]));
    j += _drho_dp[_qp] * _phi[_j][_qp] * _u_dot[_qp];
    j += _drho_dp[_qp] * _phi[_j][_qp] * _q[_qp] * _q_dot[_qp] / (_area[_qp] * _area[_qp]);
    j -= _phi[_j][_qp] * _dp_dot[_qp];
    j *= _test[_i][_qp];
  }

  return j;
}
