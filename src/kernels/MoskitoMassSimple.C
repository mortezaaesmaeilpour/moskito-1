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

#include "MoskitoMassSimple.h"

registerMooseObject("MoskitoApp", MoskitoMassSimple);

template <>
InputParameters
validParams<MoskitoMassSimple>()
{
  InputParameters params = validParams<Kernel>();

  params.addRequiredCoupledVar("flowrate", "Volumetric flowrate nonlinear variable");
  params.addRequiredCoupledVar("enthalpy", "Specific enthalpy nonlinear variable");
  params.addClassDescription("Mass conservation equation for pipe flow and "
        "it returns pressure");

  return params;
}

MoskitoMassSimple::MoskitoMassSimple(const InputParameters & parameters)
  : Kernel(parameters),
  _q_vol(coupledValue("flowrate")),
  _grad_q_vol(coupledGradient("flowrate")),
  _grad_h(coupledGradient("enthalpy")),
  _q_vol_var_number(coupled("flowrate")),
  _h_var_number(coupled("enthalpy")),
  _area(getMaterialProperty<Real>("well_area")),
  _cp(getMaterialProperty<Real>("specific_heat")),
  _well_dir(getMaterialProperty<RealVectorValue>("well_direction_vector")),
  _rho(getMaterialProperty<Real>("density")),
  _drho_dp(getMaterialProperty<Real>("drho_dp")),
  _drho_dp_2(getMaterialProperty<Real>("drho_dp_2")),
  _drho_dT(getMaterialProperty<Real>("drho_dT"))
{
}

Real
MoskitoMassSimple::computeQpResidual()
{
  RealVectorValue r = 0.0;

  r += _drho_dp[_qp] * _grad_u[_qp];
  r *= _q_vol[_qp];
  r += _grad_q_vol[_qp] * _rho[_qp];
  r *= _test[_i][_qp] / _area[_qp];

  return r * _well_dir[_qp];
}

Real
MoskitoMassSimple::computeQpJacobian()
{
  RealVectorValue j = 0.0;

  j += _drho_dp_2[_qp] * _phi[_j][_qp] * _grad_u[_qp];
  j += _drho_dp[_qp] * _grad_phi[_j][_qp];
  j *= _q_vol[_qp];
  j += _grad_q_vol[_qp] * _drho_dp[_qp] * _phi[_j][_qp];
  j *= _test[_i][_qp] / _area[_qp];

  return j * _well_dir[_qp];
}

Real
MoskitoMassSimple::computeQpOffDiagJacobian(unsigned int jvar)
{
  RealVectorValue j = 0.0;
  if (jvar == _q_vol_var_number)
  {
    j += _drho_dp[_qp] * _grad_u[_qp];
    j *= _phi[_j][_qp];
    j += _grad_phi[_j][_qp] * _rho[_qp];
    j *= _test[_i][_qp] / _area[_qp];
  }

  return j * _well_dir[_qp];
}
