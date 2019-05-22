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

#include "MoskitoHeat.h"

registerMooseObject("MoskitoApp", MoskitoHeat);

template <>
InputParameters
validParams<MoskitoHeat>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Lateral heat exchange between wellbore "
        "and formation including tubing (mandatory), insulation, liquid filled "
        "annulus and cementation");
  return params;
}

MoskitoHeat::MoskitoHeat(const InputParameters & parameters)
  : Kernel(parameters),
  _T(getMaterialProperty<Real>("temperature")),
  _rto(getMaterialProperty<Real>("radius_tubbing_outer")),
  _Uto(getMaterialProperty<Real>("thermal_resistivity_well")),
  _Twb(getMaterialProperty<Real>("temperature_formation_well_boundary"))
  {
  }

Real
MoskitoHeat::computeQpResidual()
{
  Real r = 0.0;

  r =  2.0 * PI * _rto[_qp] * _Uto[_qp];
  r *= ((_T[_qp] * gradC_to_gradR + Rankine_absol) - _Twb[_qp]);
  r /= Watt_to_Btu_per_h;
  r /= m_to_ft;
  std::cout<<"_rto = "<<_rto[_qp]<<std::endl;
  std::cout<<"_Uto = "<<_Uto[_qp]<<std::endl;
  std::cout<<"_T = "<<_T[_qp] * gradC_to_gradR + Rankine_absol<<std::endl;
  std::cout<<"_Twb = "<<_Twb[_qp]<<std::endl;
  std::cout<<"Qloss = "<<r<<std::endl;
  return - r * _test[_i][_qp];
}
