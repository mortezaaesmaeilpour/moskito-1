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

#include "MoskitoMomentum1P.h"

registerMooseObject("MoskitoApp", MoskitoMomentum1P);

template <>
InputParameters
validParams<MoskitoMomentum1P>()
{
  InputParameters params = validParams<Kernel>();

  params.addRequiredCoupledVar("pressure", "Pressure nonlinear variable");
  params.addParam<RealVectorValue>("gravity", RealVectorValue(0.0,0.0,0.0),
                                        "The gravity acceleration as a vector");
  params.addClassDescription("Momentum conservation equation for 1 phase (either"
                            " liquid or gas) pipe flow and it returns flowrate");
  return params;
}

MoskitoMomentum1P::MoskitoMomentum1P(const InputParameters & parameters)
  : Kernel(parameters),
    _grad_p(coupledGradient("pressure")),
    _p_var_number(coupled("pressure")),
    _rho(getMaterialProperty<Real>("density")),
    _drho_dp(getMaterialProperty<Real>("drho_dp")),
    _d(getMaterialProperty<Real>("well_diameter")),
    _f(getMaterialProperty<Real>("well_moody_friction")),
    _gravity(getParam<RealVectorValue>("gravity")),
    _area(getMaterialProperty<Real>("well_area")),
    _well_dir(getMaterialProperty<RealVectorValue>("well_direction_vector"))
{
}

Real
MoskitoMomentum1P::computeQpResidual()
{
  Real r = 0.0;

  r += _drho_dp[_qp] * _grad_p[_qp] * _well_dir[_qp] * _u[_qp] * _u[_qp];
  r += 2.0 * _rho[_qp] * _u[_qp] * _grad_u[_qp] * _well_dir[_qp];
  r += _f[_qp] * _rho[_qp] * _u[_qp] * fabs(_u[_qp]) / (2.0 * _d[_qp]);
  r /= (_area[_qp] * _area[_qp]);
  r += _grad_p[_qp] * _well_dir[_qp];
  r -= _rho[_qp] * _gravity * _well_dir[_qp];
  r *= _test[_i][_qp];

  return r;
}

Real
MoskitoMomentum1P::computeQpJacobian()
{
  Real j = 0.0;

  j += 2.0 * _drho_dp[_qp] * _grad_p[_qp] * _well_dir[_qp] * _phi[_j][_qp] * _u[_qp];
  j += 2.0 * _rho[_qp] * (_phi[_j][_qp]  * _grad_u[_qp] + _u[_qp] * _grad_phi[_j][_qp]) * _well_dir[_qp];
  j += 2.0 * _f[_qp] * _rho[_qp] * _phi[_j][_qp] * fabs(_u[_qp]) / (2.0 * _d[_qp]);
  j /= (_area[_qp] * _area[_qp]);
  j *= _test[_i][_qp];

  return j;
}

Real
MoskitoMomentum1P::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real j = 0.0;

  if (jvar == _p_var_number)
  {
    j += _drho_dp[_qp] * _grad_phi[_j][_qp] * _well_dir[_qp] * _u[_qp] * _u[_qp] / (_area[_qp] * _area[_qp]);
    j += _grad_phi[_j][_qp] * _well_dir[_qp];
    j *= _test[_i][_qp];
  }

  return j;
}
