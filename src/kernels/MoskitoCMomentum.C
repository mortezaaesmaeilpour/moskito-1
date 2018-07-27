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

#include "MoskitoCMomentum.h"

registerMooseObject("MoskitoApp", MoskitoCMomentum);

template <>
InputParameters
validParams<MoskitoCMomentum>()
{
  InputParameters params = validParams<Kernel>();

  params.addRequiredCoupledVar("density", "Density nonlinear variable");
  params.addCoupledVar("temperature", 273.15, "Temperature nonlinear variable");
  params.addClassDescription("Momentum equation for pipe flow based on flowrate");
  return params;
}

MoskitoCMomentum::MoskitoCMomentum(const InputParameters & parameters)
  : Kernel(parameters),
    _rho(coupledValue("density")),
    _grad_T(coupledGradient("temperature")),
    _grad_rho(coupledGradient("density")),
    _T_var_number(coupled("temperature")),
    _rho_var_number(coupled("density")),
    _d(getMaterialProperty<Real>("well_diameter")),
    _area(getMaterialProperty<Real>("well_area")),
    _dp_dT(getMaterialProperty<Real>("dp_dT")),
    _dp_drho(getMaterialProperty<Real>("dp_drho")),
    _dp_dT_2(getMaterialProperty<Real>("dp_dT_2")),
    _dp_drho_2(getMaterialProperty<Real>("dp_drho_2"))
{
}

Real
MoskitoCMomentum::computeQpResidual()
{
  Real r = 0.0;
  r += _dp_drho[_qp] * _grad_rho[_qp](0);
  r += _dp_dT[_qp] * _grad_T[_qp](0);
  r += -4.0 * _rho[_qp] * _u[_qp] * _u[_qp] / (_d[_qp] * _area[_qp] * _area[_qp]);
  r *= _test[_i][_qp];

  return r;
}

Real
MoskitoCMomentum::computeQpJacobian()
{
  Real j = 0.0;
  j += -8.0 * _rho[_qp] * _phi[_j][_qp] * _u[_qp] / (_d[_qp] * _area[_qp] * _area[_qp]) *
       _test[_i][_qp];

  return j;
}

Real
MoskitoCMomentum::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real j = 0.0;
  if (jvar == _rho_var_number)
  {
    // j += 2.0e9*(_grad_phi[_j][_qp](0) * _rho[_qp] - _phi[_j][_qp] * _grad_rho[_qp](0)) /
    // (_rho[_qp] *_rho[_qp]);

    j +=
        _dp_drho[_qp] * _grad_phi[_j][_qp](0) + _dp_drho_2[_qp] * _phi[_j][_qp] * _grad_rho[_qp](0);
    j += -4.0 * _phi[_j][_qp] * _u[_qp] * _u[_qp] / (_d[_qp] * _area[_qp] * _area[_qp]);
    j *= _test[_i][_qp];
  }

  if (jvar == _T_var_number)
  {
    j += _dp_dT[_qp] * _grad_phi[_j][_qp](0) + _dp_dT_2[_qp] * _phi[_j][_qp] * _grad_T[_qp](0);
    j *= _test[_i][_qp];
  }

  return j;
}
