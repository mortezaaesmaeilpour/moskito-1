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

#ifndef MOSKITOEOS1P_H
#define MOSKITOEOS1P_H

#include "GeneralUserObject.h"

class MoskitoEOS1P;

template <>
InputParameters validParams<MoskitoEOS1P>();

class MoskitoEOS1P : public GeneralUserObject
{
public:
  MoskitoEOS1P(const InputParameters & parameters);
  virtual ~MoskitoEOS1P();

  virtual void execute() final {}
  virtual void initialize() final {}
  virtual void finalize() final {}

  // Density from pressure and temperature (kg/m^3)
  virtual Real rho(Real pressure, Real temperature) const = 0;

  // Density from pressure and temperature and its derivatives wrt pressure and temperature
  virtual void
  drho_dpT(Real pressure, Real temperature, Real & rho, Real & drho_dp, Real & drho_dT) const = 0;

  // The second derivative of density wrt pressure, temperature and pressure-temperature
  virtual void
  drho_dpT_2(Real pressure, Real temperature, Real & drho_dp_2, Real & drho_dT_2, Real & drho_dTdp) const = 0;

  // The conversion function from temperature to specific enthalpy
  virtual Real
  T_to_h(Real temperature) const = 0;

  // The conversion function from specific enthalpy to temperature
  virtual Real
  h_to_T(Real enthalpy) const = 0;

  // specific heat at constant pressure
  Real _cp = 0;
  // thermal conductivity
  Real _lambda = 0;

protected:
  // density at reference pressure and temperature
  Real _density_ref = 0;
  // reference temperature
  Real _T_ref = 0;
  // reference pressure
  Real _P_ref = 0;
  // reference pressure
  Real _h_ref = 0;
};

#endif /* MOSKITOEOS1P_H */
