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

#ifndef MOSKITOEOS2P_H
#define MOSKITOEOS2P_H

#include "FluidProperties.h"
#include "SinglePhaseFluidProperties.h"

class MoskitoEOS2P;

template <>
InputParameters validParams<MoskitoEOS2P>();

class MoskitoEOS2P : public FluidProperties
{
public:
  MoskitoEOS2P(const InputParameters & parameters);

  virtual void GMFrac_from_p_h(const Real & p, const Real & h, Real & gmfrac, Real & T) const = 0;
  virtual Real T_from_p_h(const Real & pressure, const Real & enthalpy) const = 0;
  virtual Real h_lat(const Real & pressure, const Real & temperature) const = 0;

  const MoskitoEOS1P * liquid;
  const MoskitoEOS1P * gas;
};

#endif /* MOSKITOEOS2P_H */
