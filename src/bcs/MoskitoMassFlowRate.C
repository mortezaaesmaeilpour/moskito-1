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

#include "MoskitoMassFlowRate.h"

registerMooseObject("MoskitoApp", MoskitoMassFlowRate);

template <>
InputParameters
validParams<MoskitoMassFlowRate>()
{
  InputParameters params = validParams<NodalBC>();
  params.addRequiredParam<Real>("mass_flowrate",
        "The mass flowrate of the mixture (kg/s)");
  params.declareControllable("mass_flowrate");
  params.addRequiredParam<Real>("mixture_density",
        "The density of mixture (kg/m^3)");
  params.declareControllable("mixture_density");
  params.addClassDescription("Implements a NodalBC (Dirichlet) which calculates"
                            " mass weighted flow rate for momentum eq");
  return params;
}

MoskitoMassFlowRate::MoskitoMassFlowRate(const InputParameters & parameters)
  : NodalBC(parameters),
    _m_dot(getParam<Real>("mass_flowrate")),
    _rho_m(getParam<Real>("mixture_density"))
{
}

Real
MoskitoMassFlowRate::computeQpResidual()
{
  return _u[_qp] - (_m_dot / _rho_m);
}
