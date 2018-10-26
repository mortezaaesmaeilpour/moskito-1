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
  params.addRequiredCoupledVar("enthalpy", "Specific enthalpy nonlinear variable");
  params.addClassDescription("Momentum conservation equation for 2 phase ("
                            "liquid and gas) pipe flow and it returns flowrate");
  return params;
}

MoskitoMomentum2P::MoskitoMomentum2P(const InputParameters & parameters)
  : Kernel(parameters),
    _grad_p(coupledGradient("pressure")),
    _grad_h(coupledGradient("enthalpy")),
    _p_var_number(coupled("pressure")),
    _h_var_number(coupled("enthalpy")),
    _cp_m(getMaterialProperty<Real>("mixture_specific_heat")),
    _rho_m(getMaterialProperty<Real>("mixture_density")),
    _drho_m_dp(getMaterialProperty<Real>("drho_m_dp")),
    _drho_m_dp_2(getMaterialProperty<Real>("drho_m_dp_2")),
    _drho_m_dT(getMaterialProperty<Real>("drho_m_dT")),
    _drho_m_dT_2(getMaterialProperty<Real>("drho_m_dT_2")),
    _drho_m_dTdp(getMaterialProperty<Real>("drho_m_dTdp")),
    _drho_m_dpdT(getMaterialProperty<Real>("drho_m_dpdT")),
    _d(getMaterialProperty<Real>("well_diameter")),
    _f(getMaterialProperty<Real>("well_moody_friction")),
    _gravity(getMaterialProperty<RealVectorValue>("gravity")),
    _area(getMaterialProperty<Real>("well_area")),
    _well_dir(getMaterialProperty<RealVectorValue>("well_direction_vector")),
    _dgamma_dz(getMaterialProperty<Real>("dgamma_dz")),
    _dgamma_dz_uj_gphi(getMaterialProperty<Real>("dgamma_dz_uj_gphi")),
    _dgamma_dz_uj_phi(getMaterialProperty<Real>("dgamma_dz_uj_phi")),
    _dgamma_dz_pj_gphi(getMaterialProperty<Real>("dgamma_dz_pj_gphi")),
    _dgamma_dz_pj_phi(getMaterialProperty<Real>("dgamma_dz_pj_phi")),
    _dgamma_dz_hj_gphi(getMaterialProperty<Real>("dgamma_dz_hj_gphi")),
    _dgamma_dz_hj_phi(getMaterialProperty<Real>("dgamma_dz_hj_phi"))
{
}

Real
MoskitoMomentum2P::computeQpResidual()
{
  Real r = 0.0;

  r += _drho_m_dp[_qp] * _grad_p[_qp] * _well_dir[_qp];
  r += _drho_m_dT[_qp] * _grad_h[_qp] * _well_dir[_qp] / _cp_m[_qp];
  r *= _u[_qp] * _u[_qp];
  r += 2.0 * _rho_m[_qp] * _u[_qp] * _grad_u[_qp] * _well_dir[_qp];
  r += _f[_qp] * _rho_m[_qp] * _u[_qp] * fabs(_u[_qp]) / (2.0 * _d[_qp]);
  r /= (_area[_qp] * _area[_qp]);
  r += _dgamma_dz[_qp];
  r += _grad_p[_qp] * _well_dir[_qp];
  r -= _rho_m[_qp] * _gravity[_qp] * _well_dir[_qp];
  r *= _test[_i][_qp];

  return r;
}

Real
MoskitoMomentum2P::computeQpJacobian()
{
  Real j = 0.0;

  j += _drho_m_dp[_qp] * _grad_p[_qp] * _well_dir[_qp];
  j += _drho_m_dT[_qp] * _grad_h[_qp] * _well_dir[_qp] / _cp_m[_qp];
  j *= 2.0 * _phi[_j][_qp] * _u[_qp];
  j += 2.0 * _rho_m[_qp] * (_phi[_j][_qp]  * _grad_u[_qp] +
      _u[_qp] * _grad_phi[_j][_qp]) * _well_dir[_qp];
  j += _f[_qp] * _rho_m[_qp] * _phi[_j][_qp] * fabs(_u[_qp]) / _d[_qp];
  j /= (_area[_qp] * _area[_qp]);
  j += _dgamma_dz_uj_phi[_qp] * _phi[_j][_qp];
  j += _dgamma_dz_uj_gphi[_qp] * _grad_phi[_j][_qp] * _well_dir[_qp];
  j *= _test[_i][_qp];

  return j;
}

Real
MoskitoMomentum2P::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real j = 0.0;

  if (jvar == _p_var_number)
  {
    j += _drho_m_dp[_qp] * _grad_phi[_j][_qp] * _well_dir[_qp];
    j += _drho_m_dp_2[_qp] * _phi[_j][_qp] * _grad_p[_qp] * _well_dir[_qp];
    j += _drho_m_dTdp[_qp] * _phi[_j][_qp] * _grad_h[_qp] * _well_dir[_qp] / _cp_m[_qp];
    j *= _u[_qp] * _u[_qp];
    j += 2.0 * _drho_m_dp[_qp] * _phi[_j][_qp] * _u[_qp] * _grad_u[_qp] * _well_dir[_qp];
    j += _f[_qp] * _drho_m_dp[_qp] * _phi[_j][_qp] * _u[_qp] * fabs(_u[_qp]) / (2.0 * _d[_qp]);
    j /= (_area[_qp] * _area[_qp]);
    j += _grad_phi[_j][_qp] * _well_dir[_qp];
    j -= _drho_m_dp[_qp] * _phi[_j][_qp] * _gravity[_qp] * _well_dir[_qp];
    j += _dgamma_dz_pj_phi[_qp] * _phi[_j][_qp];
    j += _dgamma_dz_pj_gphi[_qp] * _grad_phi[_j][_qp] * _well_dir[_qp];
    j *= _test[_i][_qp];
  }

  if (jvar == _h_var_number)
  {
    j += _drho_m_dpdT[_qp] * _phi[_j][_qp] * _grad_p[_qp] * _well_dir[_qp];
    j += _drho_m_dT[_qp] * _grad_phi[_j][_qp] * _well_dir[_qp];
    j += _drho_m_dT_2[_qp] * _phi[_j][_qp] * _grad_h[_qp] * _well_dir[_qp] / _cp_m[_qp];
    j *= _u[_qp] * _u[_qp];
    j += 2.0 * _drho_m_dT[_qp] * _phi[_j][_qp] * _u[_qp] * _grad_u[_qp] * _well_dir[_qp];
    j += _f[_qp] * _drho_m_dT[_qp] * _phi[_j][_qp] * _u[_qp] * fabs(_u[_qp]) / (2.0 * _d[_qp]);
    j /= (_area[_qp] * _area[_qp]);
    j -= _drho_m_dT[_qp] * _phi[_j][_qp] * _gravity[_qp] * _well_dir[_qp];
    j /= _cp_m[_qp];
    j += _dgamma_dz_hj_phi[_qp] * _phi[_j][_qp];
    j += _dgamma_dz_hj_gphi[_qp] * _grad_phi[_j][_qp] * _well_dir[_qp];
    j *= _test[_i][_qp];
  }

  return j;
}
