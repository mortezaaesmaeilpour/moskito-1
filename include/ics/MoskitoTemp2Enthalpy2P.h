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

#ifndef MOSKITOTEMP2ENTHALPY2P_H
#define MOSKITOTEMP2ENTHALPY2P_H

// #include "InitialCondition.h"
// #include "MoskitoEOS2P.h"
#include "MoskitoWater97FluidProperties.h"
#include "InitialCondition.h"
#include "Function.h"

class MoskitoTemp2Enthalpy2P;

template <>
InputParameters validParams<MoskitoTemp2Enthalpy2P>();

class MoskitoTemp2Enthalpy2P : public InitialCondition
{
public:
  MoskitoTemp2Enthalpy2P(const InputParameters & parameters);

protected:
  virtual Real value(const Point & _point) override;

  // temperature value
  const Function & _gradT;
  // void mass fraction
  const Function & _gradvmfrac;
  // presure value
  const VariableValue & _p;

  // Userobject to equation of state
  MoskitoWater97FluidProperties * _eos_uo;
};

#endif // MOSKITOTEMPERATURETOENTHALPY1P_H
