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

#include "MoskitoFluidWell2P.h"

registerMooseObject("MoskitoApp", MoskitoFluidWell2P);

template <>
InputParameters
validParams<MoskitoFluidWell2P>()
{
  InputParameters params = validParams<MoskitoFluidWellGeneral>();
  params.addRequiredParam<UserObjectName>("eos_uo",
        "The name of the userobject for 2 phase EOS");
  return params;
}

MoskitoFluidWell2P::MoskitoFluidWell2P(const InputParameters & parameters)
  : MoskitoFluidWellGeneral(parameters),
    eos_uo(getUserObject<MoskitoEOS2P>("eos_uo")),
    _cp_m(declareProperty<Real>("mixture_specific_heat")),
    _rho_m(declareProperty<Real>("mixture_density")),
    _rho_pam(declareProperty<Real>("profile_mixture_density")),
    _drho_m_dp(declareProperty<Real>("drho_m_dp")),
    _drho_m_dp_2(declareProperty<Real>("drho_m_dp_2")),
    _drho_m_dT(declareProperty<Real>("drho_m_dT")),
    _drho_m_dT_2(declareProperty<Real>("drho_m_dT_2")),
    _drho_m_dTdp(declareProperty<Real>("drho_m_dTdp")),
    _drho_m_dpdT(declareProperty<Real>("drho_m_dpdT")),
    _alpha(declareProperty<Real>("void_fraction")),
    _u_g(declareProperty<Real>("gas velocity")),
    _u_l(declareProperty<Real>("liquid_velocity")),
    _dgamma_dz(declareProperty<Real>("dgamma_dz")),
    _dgamma_dz_uj_gphi(declareProperty<Real>("dgamma_dz_uj_gphi")),
    _dgamma_dz_uj_phi(declareProperty<Real>("dgamma_dz_uj_phi")),
    _dgamma_dz_pj_gphi(declareProperty<Real>("dgamma_dz_pj_gphi")),
    _dgamma_dz_pj_phi(declareProperty<Real>("dgamma_dz_pj_phi")),
    _dgamma_dz_hj_gphi(declareProperty<Real>("dgamma_dz_hj_gphi")),
    _dgamma_dz_hj_phi(declareProperty<Real>("dgamma_dz_hj_phi")),
    _grad_flow(coupledGradient("flowrate")),
    _grad_h(coupledGradient("enthalpy")),
    _grad_p(coupledGradient("pressure"))
{
}

void
MoskitoFluidWell2P::computeQpProperties()
{
  Real temp =0.0, c0 = 1.0, ud = 0.0;

  _alpha[_qp] = 0.0;

  eos_uo.liquid.drho_dpT(_P[_qp], _T[_qp], _rho_l, _drho_l_dp, _drho_l_dT);
  eos_uo.liquid.drho_dpT_2(_P[_qp], _T[_qp], _drho_l_dp_2, _drho_l_dT_2, temp);

  eos_uo.gas.drho_dpT(_P[_qp], _T[_qp], _rho_g, _drho_g_dp, _drho_g_dT);
  eos_uo.gas.drho_dpT_2(_P[_qp], _T[_qp], _drho_g_dp_2, _drho_g_dT_2, temp);

  _rho_m[_qp]       = _rho_g       * _alpha[_qp] + (1.0 - _alpha[_qp]) * _rho_l;
  _rho_pam[_qp]     = _rho_g * c0  * _alpha[_qp] + (1.0 - _alpha[_qp] * c0) * _rho_l;
  _drho_m_dp[_qp]   = _drho_g_dp   * _alpha[_qp] + (1.0 - _alpha[_qp]) * _drho_l_dp;
  _drho_m_dT[_qp]   = _drho_g_dT   * _alpha[_qp] + (1.0 - _alpha[_qp]) * _drho_l_dT;
  _drho_m_dp_2[_qp] = _drho_g_dp_2 * _alpha[_qp] + (1.0 - _alpha[_qp]) * _drho_l_dp_2;
  _drho_m_dT_2[_qp] = _drho_g_dT_2 * _alpha[_qp] + (1.0 - _alpha[_qp]) * _drho_l_dT_2;
  _drho_m_dpdT[_qp] = _drho_g_dpdT * _alpha[_qp] + (1.0 - _alpha[_qp]) * _drho_l_dpdT;
  _drho_m_dTdp[_qp] = _drho_g_dTdp * _alpha[_qp] + (1.0 - _alpha[_qp]) * _drho_l_dTdp;

  _u[_qp] = _flow[_qp] / _area[_qp];
  _u_g[_qp]  = c0 * _rho_m[_qp] * _u[_qp] + _rho_l * ud;
  _u_g[_qp] /= _rho_pam[_qp];
  _u_l[_qp]  = (1.0 - _alpha[_qp] * c0) * _rho_m[_qp]  * _u[_qp] - _rho_g * _alpha[_qp] * ud;
  _u_l[_qp] /= ((1.0 - _alpha[_qp]) * _rho_pam[_qp]);

  DriftFluxMomentumEq(c0, ud);

  _dia[_qp] = _d;
  _area[_qp] = PI * _d * _d / 4.0;


  // _Re[_qp] = rho_m * _dia[_qp] * fabs(_vel[_qp]) / viscosity_uo.mu(_P[_qp], _T[_qp]);

  // _lambda[_qp]  = (1.0 - (_d * _d) / std::pow(_d + _thickness , 2.0)) * _lambda0;
  // _lambda[_qp] += (_d * _d) / std::pow(_d + _thickness , 2.0) * eos_uo._lambda;

  MoskitoFluidWellGeneral::computeQpProperties();
}

void
MoskitoFluidWell2P::DriftFluxMomentumEq(Real & c0, Real & ud)
{
  /*
  All required coupled coeffient are derived below for momentum conservation
  In 2 phase momentum equation, the following rules should applied:
  1. gphi should be multiplied by (_grad_phi[j][_qp] * _well_unit_vect[_qp])
  2. phi should be multiplied by (_phi[j][_qp])
  */

  // defined and used local variables
  Real _drho_g_dz, _drho_l_dz, _drho_m_dz, _drho_pam_dz, _dc_dz;
  Real _dc_dp, _drho_pam_dp, _dc_dh, _drho_pam_dh;
  Real part1, part2, temp1, temp2;

  // parameterizing gamma equation for simpilification
  Real a, b, c, d;
  a = _alpha[_qp] / (1.0 - _alpha[_qp]);
  b = c0 - 1.0;
  c = _rho_g * _rho_l * _rho_m[_qp] / (_rho_pam[_qp] * _rho_pam[_qp]);
  d = b * _u[_qp] + ud;

  // parameterizing dgamma_dz equation for simpilification
  _drho_g_dz  = (_drho_g_dp * _grad_p[_qp] + _drho_g_dT * _grad_h[_qp] / _cp_m[_qp]) * _well_unit_vect[_qp];
  _drho_l_dz  = (_drho_l_dp * _grad_p[_qp] + _drho_l_dT * _grad_h[_qp] / _cp_m[_qp]) * _well_unit_vect[_qp];
  _drho_m_dz    =      _alpha[_qp] * _drho_g_dz + (1.0 - _alpha[_qp]     ) * _drho_l_dz;
  _drho_pam_dz  = c0 * _alpha[_qp] * _drho_g_dz + (1.0 - _alpha[_qp] * c0) * _drho_l_dz;
  part1  = _drho_g_dz * _rho_l     * _rho_m[_qp];
  part1 += _rho_g     * _drho_l_dz * _rho_m[_qp];
  part1 += _rho_g     * _rho_l     * _drho_m_dz;
  part2  = -2.0 * _drho_pam_dz * _rho_g * _rho_l * _rho_m[_qp];
  _dc_dz = (part1 * _rho_pam[_qp] + part2) / std::pow(_rho_pam[_qp] , 3.0);

  // residual for dgamma_dz in the momentum conservation
  _dgamma_dz[_qp]  = a * _dc_dz * std::pow(d, 2.0);
  _dgamma_dz[_qp] += 2.0 * a * b * c * d * _grad_flow[_qp] / _area[_qp] * _well_unit_vect[_qp];

  // diagonal jacobian of the residual wrt uj for dgamma_dz in the momentum conservation
  _dgamma_dz_uj_gphi[_qp]  = 2.0 * a * b * c * d / _area[_qp];
  _dgamma_dz_uj_phi[_qp]  = _dc_dz * d + c * b * _grad_flow[_qp] / _area[_qp] * _well_unit_vect[_qp];
  _dgamma_dz_uj_phi[_qp] *= 2.0 * a * b;

  _drho_pam_dp = c0 * _drho_g_dp * _alpha[_qp] + (1.0 - _alpha[_qp] * c0) * _drho_l_dp;
  temp1  = _drho_g_dp * _rho_l     * _rho_m[_qp];
  temp1 += _rho_g     * _drho_l_dp * _rho_m[_qp];
  temp1 += _rho_g     * _rho_l     * _drho_m_dp[_qp];
  temp2  = -2.0 * _drho_pam_dp * _rho_g * _rho_l * _rho_m[_qp];
  _dc_dp  = (temp1 * _rho_pam[_qp] + temp2) / std::pow(_rho_pam[_qp] , 3.0);

  // off diagonal jacobian of the residual wrt pj for dgamma_dz in the momentum conservation
  _dgamma_dz_pj_gphi[_qp] = a * d * d * _dc_dp;
  _dgamma_dz_pj_phi[_qp]  = _drho_g_dz * _drho_l_dp * _rho_m[_qp];
  _dgamma_dz_pj_phi[_qp] += _drho_g_dz * _rho_l     * _drho_m_dp[_qp];
  _dgamma_dz_pj_phi[_qp] += _drho_g_dp * _drho_l_dz * _rho_m[_qp];
  _dgamma_dz_pj_phi[_qp] += _rho_g     * _drho_l_dz * _drho_m_dp[_qp];
  _dgamma_dz_pj_phi[_qp] += _drho_g_dp * _rho_l     * _drho_m_dz;
  _dgamma_dz_pj_phi[_qp] += _rho_g     * _drho_l_dp * _drho_m_dz;
  _dgamma_dz_pj_phi[_qp] *= _rho_pam[_qp];
  _dgamma_dz_pj_phi[_qp] += _drho_pam_dp * part1;
  _dgamma_dz_pj_phi[_qp] -= 2.0 * _drho_m_dz * temp1;
  _dgamma_dz_pj_phi[_qp] *= _rho_pam[_qp];
  _dgamma_dz_pj_phi[_qp] -= 3.0 * _drho_pam_dp * (part1 * _rho_pam[_qp] + part2);
  _dgamma_dz_pj_phi[_qp] /= std::pow(_rho_pam[_qp] , 4.0);
  _dgamma_dz_pj_phi[_qp] *= a * d * d;
  _dgamma_dz_pj_phi[_qp] += 2.0 * a * b * _dc_dp * d * _grad_flow[_qp] / _area[_qp] * _well_unit_vect[_qp];

  _drho_pam_dh  = c0 * _drho_g_dT * _alpha[_qp] + (1.0 - _alpha[_qp] * c0) * _drho_l_dT;
  _drho_pam_dh /= _cp_m[_qp];
  temp1  = _drho_g_dT * _rho_l     * _rho_m[_qp];
  temp1 += _rho_g     * _drho_l_dT * _rho_m[_qp];
  temp1 += _rho_g     * _rho_l     * _drho_m_dT[_qp];
  temp1 /= _cp_m[_qp];
  temp2  = -2.0 * _drho_pam_dh * _rho_g * _rho_l * _rho_m[_qp];
  _dc_dh  = (temp1 * _rho_pam[_qp] + temp2) / std::pow(_rho_pam[_qp] , 3.0);

  // off diagonal jacobian of the residual wrt hj for dgamma_dz in the momentum conservation
  _dgamma_dz_hj_gphi[_qp] = a * d * d * _dc_dh;
  _dgamma_dz_hj_phi[_qp]  = _drho_g_dz * _drho_l_dT * _rho_m[_qp];
  _dgamma_dz_hj_phi[_qp] += _drho_g_dz * _rho_l     * _drho_m_dT[_qp];
  _dgamma_dz_hj_phi[_qp] += _drho_g_dT * _drho_l_dz * _rho_m[_qp];
  _dgamma_dz_hj_phi[_qp] += _rho_g     * _drho_l_dz * _drho_m_dT[_qp];
  _dgamma_dz_hj_phi[_qp] += _drho_g_dT * _rho_l     * _drho_m_dz;
  _dgamma_dz_hj_phi[_qp] += _rho_g     * _drho_l_dT * _drho_m_dz;
  _dgamma_dz_hj_phi[_qp] *= _rho_pam[_qp] / _cp_m[_qp];
  _dgamma_dz_hj_phi[_qp] += _drho_pam_dh * part1;
  _dgamma_dz_hj_phi[_qp] -= 2.0 * _drho_m_dz * temp1;
  _dgamma_dz_hj_phi[_qp] *= _rho_pam[_qp];
  _dgamma_dz_hj_phi[_qp] -= 3.0 * _drho_pam_dh * (part1 * _rho_pam[_qp] + part2);
  _dgamma_dz_hj_phi[_qp] /= std::pow(_rho_pam[_qp] , 4.0);
  _dgamma_dz_hj_phi[_qp] *= a * d * d;
  _dgamma_dz_hj_phi[_qp] += 2.0 * a * b * _dc_dh * d * _grad_flow[_qp] / _area[_qp] * _well_unit_vect[_qp];
}
