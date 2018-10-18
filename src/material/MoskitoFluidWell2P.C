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
    _drho_m_dp(declareProperty<Real>("drho_m_dp")),
    _drho_m_dp_2(declareProperty<Real>("drho_m_dp_2")),
    _drho_m_dT(declareProperty<Real>("drho_m_dT")),
    _drho_m_dT_2(declareProperty<Real>("drho_m_dT_2")),
    _drho_m_dTdp(declareProperty<Real>("drho_m_dTdp")),
    _drho_m_dpdT(declareProperty<Real>("drho_m_dpdT")),
    _alpha(declareProperty<Real>("void_fraction"))
{
}

void
MoskitoFluidWell2P::computeQpProperties()
{
  Real temp;

  eos_uo.liquid.drho_dpT(_P[_qp], _T[_qp], _rho_l, _drho_l_dp, _drho_l_dT);
  eos_uo.liquid.drho_dpT_2(_P[_qp], _T[_qp], _drho_l_dp_2, _drho_l_dT_2, temp);

  eos_uo.gas.drho_dpT(_P[_qp], _T[_qp], _rho_g, _drho_g_dp, _drho_g_dT);
  eos_uo.gas.drho_dpT_2(_P[_qp], _T[_qp], _drho_g_dp_2, _drho_g_dT_2, temp);

  _rho_m[_qp]       = _rho_g       * _alpha[_qp] + (1.0 - _alpha[_qp]) * _rho_l;
  _drho_m_dp[_qp]   = _drho_g_dp   * _alpha[_qp] + (1.0 - _alpha[_qp]) * _drho_l_dp;
  _drho_m_dT[_qp]   = _drho_g_dT   * _alpha[_qp] + (1.0 - _alpha[_qp]) * _drho_l_dT;
  _drho_m_dp_2[_qp] = _drho_g_dp_2 * _alpha[_qp] + (1.0 - _alpha[_qp]) * _drho_l_dp_2;
  _drho_m_dT_2[_qp] = _drho_g_dT_2 * _alpha[_qp] + (1.0 - _alpha[_qp]) * _drho_l_dT_2;
  _drho_m_dpdT[_qp] = _drho_g_dpdT * _alpha[_qp] + (1.0 - _alpha[_qp]) * _drho_l_dpdT;
  _drho_m_dTdp[_qp] = _drho_g_dTdp * _alpha[_qp] + (1.0 - _alpha[_qp]) * _drho_l_dTdp;

  _dia[_qp] = _d;
  _area[_qp] = PI * _d * _d / 4.0;

  // _vel[_qp] = _flow[_qp] / _area[_qp];
  // _Re[_qp] = rho_m * _dia[_qp] * fabs(_vel[_qp]) / viscosity_uo.mu(_P[_qp], _T[_qp]);

  // _lambda[_qp]  = (1.0 - (_d * _d) / std::pow(_d + _thickness , 2.0)) * _lambda0;
  // _lambda[_qp] += (_d * _d) / std::pow(_d + _thickness , 2.0) * eos_uo._lambda;

  MoskitoFluidWellGeneral::computeQpProperties();
}
