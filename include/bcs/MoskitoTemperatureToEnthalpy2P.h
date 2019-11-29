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

#ifndef MOSKITOTEMPERATURETOENTHALPY2P_H
#define MOSKITOTEMPERATURETOENTHALPY2P_H

#include "NodalBC.h"
#include "MoskitoEOS2P.h"
#include "MoskitoWater97FluidProperties.h"

class MoskitoTemperatureToEnthalpy2P;

template <>
InputParameters validParams<MoskitoTemperatureToEnthalpy2P>();

class MoskitoTemperatureToEnthalpy2P : public NodalBC
{
public:
  MoskitoTemperatureToEnthalpy2P(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  // temperature value
  const Real & _T;
  // void mass fraction
  const Real & _vmfrac;
  // presure value
  const VariableValue & _p;

  // Userobject to equation of state
  MoskitoWater97FluidProperties * _eos_uo;
};

#endif // MOSKITOTEMPERATURETOENTHALPY1P_H
