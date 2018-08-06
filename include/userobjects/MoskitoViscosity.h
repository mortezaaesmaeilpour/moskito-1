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

#ifndef MOSKITOVISCOSITY_H
#define MOSKITOVISCOSITY_H

#include "GeneralUserObject.h"

class MoskitoViscosity;

template <>
InputParameters validParams<MoskitoViscosity>();

class MoskitoViscosity : public GeneralUserObject
{
public:
  MoskitoViscosity(const InputParameters & parameters);
  virtual ~MoskitoViscosity();

  virtual void execute() final {}
  virtual void initialize() final {}
  virtual void finalize() final {}

  // viscosity from pressure and temperature (pa.s)
  virtual Real mu(Real pressure, Real temperature) const = 0;

  // viscosity from pressure and temperature and its derivatives wrt pressure and temperature
  virtual void
  dmu_dpT(Real pressure, Real temperature, Real & rho, Real & drho_dp, Real & drho_dT) const = 0;
};

#endif /* MOSKITOVISCOSITY_H */
