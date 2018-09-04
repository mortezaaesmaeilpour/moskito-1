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

#include "MoskitoMomentum2P.h"

registerMooseObject("MoskitoApp", MoskitoMomentum2P);

template <>
InputParameters
validParams<MoskitoMomentum2P>()
{
  InputParameters params = validParams<Kernel>();

  params.addRequiredCoupledVar("pressure", "Pressure nonlinear variable");
  params.addRequiredCoupledVar("void_fraction", "Gas saturation nonlinear variable");
  params.addClassDescription("Momentum equation for pipe flow based on flowrate");
  params.addParam<RealVectorValue>("gravity", RealVectorValue(0.0,0.0,0.0), "The gravity acceleration as a vector");
  return params;
}

MoskitoMomentum2P::MoskitoMomentum2P(const InputParameters & parameters)
  : Kernel(parameters),
    _alpha(coupledValue("void_fraction")),
    _grad_p(coupledGradient("pressure")),
    _grad_alpha(coupledGradient("void_fraction")),
    _p_var_number(coupled("pressure")),
    _alpha_var_number(coupled("void_fraction")),
    _rho_g(getMaterialProperty<Real>("gas_density")),
    _rho_l(getMaterialProperty<Real>("liquid_density")),
    _drho_g_dp(getMaterialProperty<Real>("drho_g_dp")),
    _drho_l_dp(getMaterialProperty<Real>("drho_l_dp")),
    _d(getMaterialProperty<Real>("well_diameter")),
    _f(getMaterialProperty<Real>("well_moody_friction")),
    _gravity(getParam<RealVectorValue>("gravity")),
    _area(getMaterialProperty<Real>("well_area")),
    _well_dir(getMaterialProperty<RealVectorValue>("well_direction_vector"))
{
}

Real
MoskitoMomentum2P::computeQpResidual()
{
  Real r = 0.0, rho_m = 0.0, drho_m_dl = 0.0;

  rho_m = _alpha[_qp] * _rho_g[_qp] + (1.0 - _alpha[_qp]) * _rho_l[_qp];

  drho_m_dl += _grad_alpha[_qp] * _well_dir[_qp] * (_rho_g[_qp] - _rho_l[_qp]);
  drho_m_dl += _grad_p[_qp] * _well_dir[_qp] * ( _alpha[_qp] * (_drho_g_dp[_qp]
                                  - _drho_l_dp[_qp]) + _drho_l_dp[_qp]);

  r += drho_m_dl * _u[_qp] * _u[_qp] / (_area[_qp] * _area[_qp]);
  r += 2.0 * rho_m * _u[_qp] * _grad_u[_qp]  * _well_dir[_qp] / (_area[_qp] * _area[_qp]);
  r += _grad_p[_qp] * _well_dir[_qp];
  r += _f[_qp] * rho_m * _u[_qp] * fabs(_u[_qp]) / (2.0 * _d[_qp] * _area[_qp] * _area[_qp]);
  r -= rho_m * _gravity * _well_dir[_qp];

  return r * _test[_i][_qp];
}

Real
MoskitoMomentum2P::computeQpJacobian()
{
  Real j = 0.0, rho_m = 0.0, drho_m_dl = 0.0;

  rho_m = _alpha[_qp] * _rho_g[_qp] + (1.0 - _alpha[_qp]) * _rho_l[_qp];

  drho_m_dl += _grad_alpha[_qp] * _well_dir[_qp] * (_rho_g[_qp] - _rho_l[_qp]);
  drho_m_dl += _grad_p[_qp] * _well_dir[_qp] * ( _alpha[_qp] * (_drho_g_dp[_qp]
                                  - _drho_l_dp[_qp]) + _drho_l_dp[_qp]);

  j += 2.0 * drho_m_dl / (_area[_qp] * _area[_qp]) * _phi[_j][_qp] * _u[_qp];
  j += 2.0 * rho_m / (_area[_qp] * _area[_qp]) * _well_dir[_qp] * (_phi[_j][_qp]
                                * _grad_u[_qp] + _u[_qp] * _grad_phi[_j][_qp]);
  j += 2.0 * _f[_qp] * rho_m * _phi[_j][_qp] * fabs(_u[_qp]) /
       (2.0 * _d[_qp] * _area[_qp] * _area[_qp]);

  return j * _test[_i][_qp];
}

Real
MoskitoMomentum2P::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real j = 0.0, rho_m = 0.0, drho_m_dl = 0.0;

  if (jvar == _p_var_number)
  {
    drho_m_dl += _grad_phi[_j][_qp] * _well_dir[_qp] * ( _alpha[_qp] * (_drho_g_dp[_qp]
                                    - _drho_l_dp[_qp]) + _drho_l_dp[_qp]);

    j += drho_m_dl * _u[_qp] * _u[_qp] / (_area[_qp] * _area[_qp]);
    j += _grad_phi[_j][_qp] * _well_dir[_qp];
  }

  if (jvar == _alpha_var_number)
  {
    rho_m = _phi[_j][_qp] * _rho_g[_qp] + (1.0 - _phi[_j][_qp]) * _rho_l[_qp];

    drho_m_dl += _grad_phi[_j][_qp] * _well_dir[_qp] * (_rho_g[_qp] - _rho_l[_qp]);
    drho_m_dl += _grad_p[_qp] * _well_dir[_qp] * ( _phi[_j][_qp] * (_drho_g_dp[_qp]
                                    - _drho_l_dp[_qp]));

    j += drho_m_dl * _u[_qp] * _u[_qp] / (_area[_qp] * _area[_qp]);
    j += 2.0 * rho_m * _u[_qp] * _grad_u[_qp] * _well_dir[_qp] / (_area[_qp] * _area[_qp]);
    j += _f[_qp] * rho_m * _u[_qp] * fabs(_u[_qp]) / (2.0 * _d[_qp] * _area[_qp] * _area[_qp]);
    j -= rho_m * _gravity * _well_dir[_qp];
  }

  return j * _test[_i][_qp];
}
