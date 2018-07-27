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

#ifndef MOSKITOEOS_H
#define MOSKITOEOS_H

#include "GeneralUserObject.h"

class MoskitoEOS;

template <>
InputParameters validParams<MoskitoEOS>();

class MoskitoEOS : public GeneralUserObject
{
public:
  MoskitoEOS(const InputParameters & parameters);
  virtual ~MoskitoEOS();

  virtual void execute() final {}
  virtual void initialize() final {}
  virtual void finalize() final {}

  /// Pressure from density and temperature (Pa)
  virtual Real p(Real density, Real temperature) const = 0;

  /// Pressure from density and temperature and its derivatives wrt density and temperature
  virtual void
  dp_drhoT(Real density, Real temperature, Real & pressure, Real & dp_drho, Real & dp_dT) const = 0;

  virtual void
  dp_drhoT_2(Real density, Real temperature, Real & dp_drho_2, Real & dp_dT_2) const = 0;
};

#endif /* MOSKITOEOS_H */
