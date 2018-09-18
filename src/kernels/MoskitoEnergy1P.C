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

#include "MoskitoEnergy1P.h"

registerMooseObject("MoskitoApp", MoskitoEnergy1P);

template <>
InputParameters
validParams<MoskitoEnergy1P>()
{
  InputParameters params = validParams<Kernel>();

  params.addRequiredCoupledVar("flowrate", "Volumetric flowrate nonlinear variable");
  params.addRequiredCoupledVar("pressure", "Pressure nonlinear variable");
  params.addClassDescription("Energy conservation equation for 1 phase (either liquid or"
                                      " gas) pipe flow and it returns enthalpy");

  return params;
}

MoskitoEnergy1P::MoskitoEnergy1P(const InputParameters & parameters)
  : Kernel(parameters),
  _q_vol(coupledValue("flowrate")),
  _grad_q_vol(coupledGradient("flowrate")),
  _grad_p(coupledGradient("pressure")),
  _q_vol_var_number(coupled("flowrate")),
  _area(getMaterialProperty<Real>("well_area")),
  _well_dir(getMaterialProperty<RealVectorValue>("well_direction_vector")),
  _cp(getMaterialProperty<Real>("specific_heat")),
  _rho(getMaterialProperty<Real>("density")),
  _drho_dp(getMaterialProperty<Real>("drho_dp")),
  _drho_dp_2(getMaterialProperty<Real>("drho_dp_2")),
  _drho_dT(getMaterialProperty<Real>("drho_dT")),
  _drho_dT_2(getMaterialProperty<Real>("drho_dT_2")),
  _drho_dTdp(getMaterialProperty<Real>("drho_dTdp")),
  _gravity(getMaterialProperty<RealVectorValue>("gravity"))
{
}

Real
MoskitoEnergy1P::computeQpResidual()
{
  Real r = 0.0, V = 0.0, grad_V = 0.0, grad_rho_V = 0.0;

  V = _q_vol[_qp] / _area[_qp];

  grad_V = _grad_q_vol[_qp] * _well_dir[_qp] / _area[_qp];

  grad_rho_V += _drho_dp[_qp] * _grad_p[_qp] * _well_dir[_qp];
  grad_rho_V += _drho_dT[_qp] * _grad_u[_qp] * _well_dir[_qp] / _cp[_qp];
  grad_rho_V *= V;
  grad_rho_V +=  grad_V * _rho[_qp];

  r += grad_rho_V * (_u[_qp] + V * V / 2.0);
  r += _rho[_qp] * V * _grad_u[_qp] * _well_dir[_qp];
  r += _rho[_qp] * V * V * grad_V;
  r -= _rho[_qp] * V * _gravity[_qp] * _well_dir[_qp];


  return r * _test[_i][_qp];
}

Real
MoskitoEnergy1P::computeQpJacobian()
{
  Real j = 0.0, V = 0.0, grad_V = 0.0, grad_rho_V = 0.0, grad_rho_V_Uj = 0.0;

  V = _q_vol[_qp] / _area[_qp];

  grad_V = _grad_q_vol[_qp] * _well_dir[_qp] / _area[_qp];

  grad_rho_V += _drho_dp[_qp] * _grad_p[_qp] * _well_dir[_qp];
  grad_rho_V += _drho_dT[_qp] * _grad_u[_qp] * _well_dir[_qp] / _cp[_qp];
  grad_rho_V *= V;
  grad_rho_V +=  grad_V * _rho[_qp];

  grad_rho_V_Uj += _drho_dTdp[_qp] * _phi[_j][_qp] * _grad_p[_qp] * _well_dir[_qp] / _cp[_qp];
  grad_rho_V_Uj += _drho_dT_2[_qp] * _phi[_j][_qp] * _grad_u[_qp] * _well_dir[_qp] / (_cp[_qp] * _cp[_qp]);
  grad_rho_V_Uj += _drho_dT[_qp] * _grad_phi[_j][_qp] * _well_dir[_qp] / _cp[_qp];
  grad_rho_V_Uj *= V;
  grad_rho_V_Uj +=  grad_V * _drho_dT[_qp] * _phi[_j][_qp] / _cp[_qp];

  j += grad_rho_V_Uj * (_u[_qp] + V * V / 2.0);
  j += grad_rho_V * _phi[_j][_qp];
  j += _rho[_qp] * V * _grad_phi[_j][_qp] * _well_dir[_qp];
  j += _drho_dT[_qp] * _phi[_j][_qp] / _cp[_qp] * V * _grad_u[_qp] * _well_dir[_qp];
  j += _drho_dT[_qp] * _phi[_j][_qp] / _cp[_qp] * V * V * grad_V;
  j -= _drho_dT[_qp] * _phi[_j][_qp] / _cp[_qp] * V * _gravity[_qp] * _well_dir[_qp];

  return j * _test[_i][_qp];
}

Real
MoskitoEnergy1P::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real j = 0.0, V = 0.0, grad_V = 0.0, grad_rho_V = 0.0;

  if (jvar == _q_vol_var_number)
  {
    Real V_Vj = 0.0, grad_V_Vj = 0.0, grad_rho_V_Vj = 0.0;

    V_Vj = _phi[_j][_qp] / _area[_qp];

    grad_V_Vj = _grad_phi[_j][_qp] * _well_dir[_qp] / _area[_qp];

    grad_rho_V_Vj += _drho_dp[_qp] * _grad_p[_qp] * _well_dir[_qp];
    grad_rho_V_Vj += _drho_dT[_qp] * _grad_u[_qp] * _well_dir[_qp] / _cp[_qp];
    grad_rho_V_Vj *= V_Vj;
    grad_rho_V_Vj +=  grad_V_Vj * _rho[_qp];

    V = _q_vol[_qp] / _area[_qp];

    grad_V = _grad_q_vol[_qp] * _well_dir[_qp] / _area[_qp];

    grad_rho_V += _drho_dp[_qp] * _grad_p[_qp] * _well_dir[_qp];
    grad_rho_V += _drho_dT[_qp] * _grad_u[_qp] * _well_dir[_qp] / _cp[_qp];
    grad_rho_V *= V;
    grad_rho_V +=  grad_V * _rho[_qp];

    j += grad_rho_V_Vj * (_u[_qp] + V * V / 2.0);
    j += grad_rho_V * V * V_Vj;
    j += _rho[_qp] * V_Vj * _grad_u[_qp] * _well_dir[_qp];
    j += _rho[_qp] * V * V * grad_V_Vj;
    j += 2.0 * _rho[_qp] * V * V_Vj * grad_V;
    j -= _rho[_qp] * V_Vj * _gravity[_qp] * _well_dir[_qp];
  }

  if (jvar == _p_var_number)
  {
    Real grad_rho_V_Pj = 0.0;

    V = _q_vol[_qp] / _area[_qp];

    grad_V = _grad_q_vol[_qp] * _well_dir[_qp] / _area[_qp];

    grad_rho_V_Pj += _drho_dp_2[_qp] * _phi[_j][_qp] * _grad_p[_qp] * _well_dir[_qp];
    grad_rho_V_Pj += _drho_dp[_qp] * _grad_phi[_j][_qp] * _well_dir[_qp];
    grad_rho_V_Pj *= V;
    grad_rho_V_Pj +=  grad_V * _drho_dp[_qp] * _phi[_j][_qp];

    j += grad_rho_V_Pj * (_u[_qp] + V * V / 2.0);
    j += _drho_dp[_qp] * _phi[_j][_qp] * V * _grad_u[_qp] * _well_dir[_qp];
    j += _drho_dp[_qp] * _phi[_j][_qp] * V * V * grad_V;
    j -= _drho_dp[_qp] * _phi[_j][_qp] * V * _gravity[_qp] * _well_dir[_qp];
  }

  return j * _test[_i][_qp];
}
