/**************************************************************************/
/*  TIGER - Hydro-thermal sImulator GEothermal Reservoirs                 */
/*                                                                        */
/*  Copyright (C) 2017 by Maziar Gholami Korzani                          */
/*  Karlsruhe Institute of Technology, Institute of Applied Geosciences   */
/*  Division of Geothermal Research                                       */
/*                                                                        */
/*  This file is part of TIGER App                                        */
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

#include "MoskitoCMomentumMaterial.h"

registerMooseObject("MoskitoApp", MoskitoCMomentumMaterial);

template <>
InputParameters
validParams<MoskitoCMomentumMaterial>()
{
  InputParameters params = validParams<Material>();

  params.addRequiredCoupledVar("density", "Density nonlinear variable (kg/m^3)");
  params.addCoupledVar("temperature", 273.15, "Temperature nonlinear variable (K)");

  params.addRequiredParam<Real>("well_diameter", "Well diameter (m)");
  params.addRequiredParam<UserObjectName>("eos_UO", "The name of the userobject for EOS");

  return params;
}

MoskitoCMomentumMaterial::MoskitoCMomentumMaterial(const InputParameters & parameters)
  : Material(parameters),
    _dia(declareProperty<Real>("well_diameter")),
    _area(declareProperty<Real>("well_area")),
    _p(declareProperty<Real>("pressure")),
    _dp_drho(declareProperty<Real>("dp_drho")),
    _dp_dT(declareProperty<Real>("dp_dT")),
    _dp_drho_2(declareProperty<Real>("dp_drho_2")),
    _dp_dT_2(declareProperty<Real>("dp_dT_2")),
    _eos_UO(getUserObject<MoskitoEOS>("eos_UO")),
    _rho(coupledValue("density")),
    _T(coupledValue("temperature")),
    _d(getParam<Real>("well_diameter"))
{
}

void
MoskitoCMomentumMaterial::computeQpProperties()
{
  _dia[_qp] = _d;
  _area[_qp] = PI * _d * _d / 4.0;

  Real p, dp_drho, dp_dT;

  _eos_UO.dp_drhoT(_rho[_qp], _T[_qp], p, dp_drho, dp_dT);
  _p[_qp] = p;
  _dp_drho[_qp] = dp_drho;
  _dp_dT[_qp] = dp_dT;

  _eos_UO.dp_drhoT_2(_rho[_qp], _T[_qp], dp_drho, dp_dT);
  _dp_drho_2[_qp] = dp_drho;
  _dp_dT_2[_qp] = dp_dT;

}
