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

#include "MoskitoMass2P.h"

registerMooseObject("MoskitoApp", MoskitoMass2P);

template <>
InputParameters
validParams<MoskitoMass2P>()
{
  InputParameters params = validParams<Kernel>();

  params.addRequiredCoupledVar("flow_rate", "Volumetric flow rate nonlinear variable");
  params.addRequiredCoupledVar("void_fraction", "Gas saturation nonlinear variable");

  return params;
}

MoskitoMass2P::MoskitoMass2P(const InputParameters & parameters)
  : Kernel(parameters),
  _q_vol(coupledValue("flow_rate")),
  _alpha(coupledValue("void_fraction")),
  _grad_q_vol(coupledGradient("flow_rate")),
  _grad_alpha(coupledGradient("void_fraction")),
  _q_vol_var_number(coupled("flow_rate")),
  _alpha_var_number(coupled("void_fraction")),
  _area(getMaterialProperty<Real>("well_area")),
  _well_dir(getMaterialProperty<RealVectorValue>("well_direction_vector")),
  _rho_g(getMaterialProperty<Real>("gas_density")),
  _rho_l(getMaterialProperty<Real>("liquid_density")),
  _drho_g_dp(getMaterialProperty<Real>("drho_g_dp")),
  _drho_l_dp(getMaterialProperty<Real>("drho_l_dp")),
  _drho_g_dp_2(getMaterialProperty<Real>("drho_g_dp_2")),
  _drho_l_dp_2(getMaterialProperty<Real>("drho_l_dp_2"))
{
}

Real
MoskitoMass2P::computeQpResidual()
{
  // grad_rho * v + rho(p) * grad_v
  // drho_dp * grad_p * v + rho(p) * grad_v
  Real r = 0.0, rho_m = 0.0;

  rho_m = _alpha[_qp] * _rho_g[_qp] + (1.0 - _alpha[_qp]) * _rho_l[_qp];
  r += _grad_alpha[_qp] * _well_dir[_qp] * (_rho_g[_qp] - _rho_l[_qp]);
  r += _grad_u[_qp] * _well_dir[_qp] * ( _alpha[_qp] * (_drho_g_dp[_qp]
                                  - _drho_l_dp[_qp]) + _drho_l_dp[_qp]);
  r *= _q_vol[_qp];
  r += _grad_q_vol[_qp] * _well_dir[_qp] * rho_m;

  return r / _area[_qp] * _test[_i][_qp];
}

Real
MoskitoMass2P::computeQpJacobian()
{
  Real j = 0.0, rho_m_uj = 0.0;

  rho_m_uj  = _alpha[_qp] * _drho_g_dp[_qp] * _phi[_j][_qp];
  rho_m_uj += (1.0 - _alpha[_qp]) * _drho_l_dp[_qp] * _phi[_j][_qp];
  j += _grad_alpha[_qp] * _well_dir[_qp] * (_drho_g_dp[_qp] - _drho_l_dp[_qp]) * _phi[_j][_qp];
  j += _grad_phi[_j][_qp] * _well_dir[_qp] * ( _alpha[_qp] * (_drho_g_dp[_qp]
                                          - _drho_l_dp[_qp]) + _drho_l_dp[_qp]);
  j += _grad_u[_qp] * _well_dir[_qp] * ( _alpha[_qp] * (_drho_g_dp_2[_qp]
                      - _drho_l_dp_2[_qp]) + _drho_l_dp_2[_qp]) * _phi[_j][_qp];
  j *= _q_vol[_qp];
  j += _grad_q_vol[_qp] * _well_dir[_qp] * rho_m_uj;

  return j / _area[_qp] * _test[_i][_qp];
}

Real
MoskitoMass2P::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real j = 0.0, rho_m = 0.0;
  if (jvar == _q_vol_var_number)
  {
    rho_m = _alpha[_qp] * _rho_g[_qp] + (1.0 - _alpha[_qp]) * _rho_l[_qp];
    j += _grad_alpha[_qp] * _well_dir[_qp] * (_rho_g[_qp] - _rho_l[_qp]);
    j += _grad_u[_qp] * _well_dir[_qp] * ( _alpha[_qp] * (_drho_g_dp[_qp]
                                    - _drho_l_dp[_qp]) + _drho_l_dp[_qp]);
    j *= _phi[_j][_qp];
    j += _grad_phi[_j][_qp] * _well_dir[_qp] * rho_m;
  }

  if (jvar == _alpha_var_number)
  {
    rho_m = _phi[_j][_qp] * _rho_g[_qp] - _phi[_j][_qp] * _rho_l[_qp];
    j += _grad_phi[_j][_qp] * _well_dir[_qp] * (_rho_g[_qp] - _rho_l[_qp]);
    j += _grad_u[_qp] * _well_dir[_qp] * ( _phi[_j][_qp] * (_drho_g_dp[_qp]
                                    - _drho_l_dp[_qp]));
    j *= _q_vol[_qp];
    j += _grad_q_vol[_qp] * _well_dir[_qp] * rho_m;
  }

  return j / _area[_qp] * _test[_i][_qp];
}
