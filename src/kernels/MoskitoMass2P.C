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

  params.addRequiredCoupledVar("flowrate", "Volumetric flowrate nonlinear variable");
  params.addRequiredCoupledVar("enthalpy", "Specific enthalpy nonlinear variable");
  params.addClassDescription("Mass conservation equation for 2 phase (liquid and"
                                      " gas) pipe flow and it returns pressure");

  return params;
}

MoskitoMass2P::MoskitoMass2P(const InputParameters & parameters)
  : Kernel(parameters),
  _q_vol(coupledValue("flowrate")),
  _grad_q_vol(coupledGradient("flowrate")),
  _grad_h(coupledGradient("enthalpy")),
  _q_vol_var_number(coupled("flowrate")),
  _h_var_number(coupled("enthalpy")),
  _area(getMaterialProperty<Real>("well_area")),
  _cp(getMaterialProperty<Real>("specific_heat")),
  _well_dir(getMaterialProperty<RealVectorValue>("well_direction_vector")),
  _rho_m(getMaterialProperty<Real>("mixture_density")),
  _drho_m_dp(getMaterialProperty<Real>("drho_m_dp")),
  _drho_m_dp_2(getMaterialProperty<Real>("drho_m_dp_2")),
  _drho_m_dT(getMaterialProperty<Real>("drho_m_dT")),
  _drho_m_dT_2(getMaterialProperty<Real>("drho_m_dT_2")),
  _drho_m_dTdp(getMaterialProperty<Real>("drho_m_dTdp")),
  _drho_m_dpdT(getMaterialProperty<Real>("drho_m_dpdT"))
{
}

Real
MoskitoMass2P::computeQpResidual()
{
  Real r = 0.0;

  r += _drho_m_dp[_qp] * _grad_u[_qp] * _well_dir[_qp];
  r += _drho_m_dT[_qp] * _grad_h[_qp] * _well_dir[_qp] / _cp[_qp];
  r *= _q_vol[_qp];
  r += _grad_q_vol[_qp] * _well_dir[_qp] * _rho_m[_qp];

  return r * _test[_i][_qp] / _area[_qp];
}

Real
MoskitoMass2P::computeQpJacobian()
{
  Real j = 0.0;

  j += _drho_m_dp_2[_qp] * _phi[_j][_qp] * _grad_u[_qp] * _well_dir[_qp];
  j += _drho_m_dp[_qp] * _grad_phi[_j][_qp] * _well_dir[_qp];
  j += _drho_m_dTdp[_qp] * _phi[_j][_qp] * _grad_h[_qp] * _well_dir[_qp] / _cp[_qp];
  j *= _q_vol[_qp];
  j += _grad_q_vol[_qp] * _drho_m_dp[_qp] * _phi[_j][_qp] * _well_dir[_qp];
  j *= _test[_i][_qp] / _area[_qp];

  return j;
}

Real
MoskitoMass2P::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real j = 0.0;
  if (jvar == _q_vol_var_number)
  {
    j += _drho_m_dp[_qp] * _grad_u[_qp] * _well_dir[_qp];
    j += _drho_m_dT[_qp] * _grad_h[_qp] * _well_dir[_qp] / _cp[_qp];
    j *= _phi[_j][_qp];
    j += _grad_phi[_j][_qp] * _rho_m[_qp] * _well_dir[_qp];
    j *= _test[_i][_qp] / _area[_qp];
  }

  if (jvar == _h_var_number)
  {
    j += _drho_m_dpdT[_qp] * _phi[_j][_qp] * _grad_u[_qp] * _well_dir[_qp];
    j += _drho_m_dT_2[_qp] * _phi[_j][_qp] * _grad_h[_qp] * _well_dir[_qp] / _cp[_qp];
    j += _drho_m_dT[_qp] * _grad_phi[_j][_qp] * _well_dir[_qp];
    j *= _q_vol[_qp];
    j += _grad_q_vol[_qp] * _well_dir[_qp] * _drho_m_dT[_qp] * _phi[_j][_qp];
    j *= _test[_i][_qp] / (_area[_qp] * _cp[_qp]);
  }

  return j;
}
