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

  params.addCoupledVar("void_fraction", 0.0,"Void fraction nonlinear variable (-)");

  return params;
}

MoskitoFluidWell2P::MoskitoFluidWell2P(const InputParameters & parameters)
  : MoskitoFluidWellGeneral(parameters),
    _rho_g(declareProperty<Real>("gas_density")),
    _rho_l(declareProperty<Real>("liquid_density")),
    _drho_g_dp(declareProperty<Real>("drho_g_dp")),
    _drho_l_dp(declareProperty<Real>("drho_l_dp")),
    _drho_g_dp_2(declareProperty<Real>("drho_g_dp_2")),
    _drho_l_dp_2(declareProperty<Real>("drho_l_dp_2")),
    _alpha(coupledValue("void_fraction"))
{
}

void
MoskitoFluidWell2P::computeQpProperties()
{
  Real rho_m, temp, temp2;

  _eos_UO.drho_dpT(_P[_qp], _T[_qp], _rho_l[_qp], _drho_l_dp[_qp], temp);
  _eos_UO.drho_dpT_2(_P[_qp], _T[_qp], _drho_l_dp_2[_qp], temp, temp2);

  _rho_g[_qp] = 0.0;
  _drho_g_dp[_qp] = 0.0;
  _drho_g_dp_2[_qp] = 0.0;

  rho_m = _rho_g[_qp] * _alpha[_qp] + (1.0 - _alpha[_qp]) * _rho_l[_qp];

  _dia[_qp] = _d;
  _area[_qp] = PI * _d * _d / 4.0;

  _vel[_qp] = _flow[_qp] / _area[_qp];
  _Re[_qp] = rho_m * _dia[_qp] * fabs(_vel[_qp]) / _viscosity_UO.mu(_P[_qp], _T[_qp]);

  MoskitoFluidWellGeneral::computeQpProperties();
}
