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

#include "MoskitoMass1P.h"

registerMooseObject("MoskitoApp", MoskitoMass1P);

template <>
InputParameters
validParams<MoskitoMass1P>()
{
  InputParameters params = validParams<Kernel>();

  params.addRequiredCoupledVar("flow_rate", "Volumetric flow rate nonlinear variable");
  params.addClassDescription("Mass conservation equation for 1 phase (either liquid or"
                                      " gas) pipe flow and it returns pressure");

  return params;
}

MoskitoMass1P::MoskitoMass1P(const InputParameters & parameters)
  : Kernel(parameters),
  _q_vol(coupledValue("flow_rate")),
  _grad_q_vol(coupledGradient("flow_rate")),
  _q_vol_var_number(coupled("flow_rate")),
  _area(getMaterialProperty<Real>("well_area")),
  _well_dir(getMaterialProperty<RealVectorValue>("well_direction_vector")),
  _rho(getMaterialProperty<Real>("density")),
  _drho_dp(getMaterialProperty<Real>("drho_dp")),
  _drho_dp_2(getMaterialProperty<Real>("drho_dp_2"))
{
}

Real
MoskitoMass1P::computeQpResidual()
{
  RealVectorValue r = 0.0;

  r += _drho_dp[_qp] * _grad_u[_qp] * _q_vol[_qp];
  r += _grad_q_vol[_qp] * _rho[_qp];

  return r * _test[_i][_qp] * _well_dir[_qp] / _area[_qp];
}

Real
MoskitoMass1P::computeQpJacobian()
{
  RealVectorValue j = 0.0;

  j += _drho_dp_2[_qp] * _phi[_j][_qp] * _grad_u[_qp];
  j += _drho_dp[_qp] * _grad_phi[_j][_qp];
  j *= _q_vol[_qp];
  j += _grad_q_vol[_qp] * _drho_dp[_qp] * _phi[_j][_qp];

  return j * _test[_i][_qp] * _well_dir[_qp] / _area[_qp];
}

Real
MoskitoMass1P::computeQpOffDiagJacobian(unsigned int jvar)
{
  RealVectorValue j = 0.0;
  if (jvar == _q_vol_var_number)
  {
    j += _drho_dp[_qp] * _grad_u[_qp] * _phi[_j][_qp];
    j += _grad_phi[_j][_qp] * _rho[_qp];
  }

  return j * _test[_i][_qp] * _well_dir[_qp] / _area[_qp];
}
