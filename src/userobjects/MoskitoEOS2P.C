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

#include "MoskitoEOS2P.h"

template <>
InputParameters
validParams<MoskitoEOS2P>()
{
  InputParameters params = validParams<GeneralUserObject>();

  return params;
}

MoskitoEOS2P::MoskitoEOS2P(const InputParameters & parameters)
  : FluidProperties(parameters)
{
}

Real
MoskitoEOS2P::rho_m_from_p_T(
      const Real & pressure, const Real & temperature, const Real & vmfrac, const unsigned int & phase) const
{
  Real rho = 0.0;

  switch (phase)
  {
    case 0:
      rho  = rho_l_from_p_T(pressure, temperature, phase);
      break;

    case 1:
      rho  = rho_g_from_p_T(pressure, temperature, phase);
      break;

    case 2:
      rho  = rho_g_from_p_T(pressure, temperature, phase) * vmfrac;
      rho += rho_l_from_p_T(pressure, temperature, phase) * (1.0 - vmfrac);
      break;

    case 3:
      rho  = rho_l_from_p_T(pressure, temperature, phase);
      break;
  }

  return rho;
}
