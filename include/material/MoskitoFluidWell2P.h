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

#ifndef MOSKITOFLUIDWELL2P_H
#define MOSKITOFLUIDWELL2P_H

#include "MoskitoFluidWellGeneral.h"

class MoskitoFluidWell2P;

template <>
InputParameters validParams<MoskitoFluidWell2P>();

class MoskitoFluidWell2P : public MoskitoFluidWellGeneral
{
public:
  MoskitoFluidWell2P(const InputParameters & parameters);
  virtual void computeQpProperties() override;

protected:
  // density of gas
  MaterialProperty<Real> & _rho_g;
  // density of liquid
  MaterialProperty<Real> & _rho_l;
  // The first derivative of gas density wrt pressure
  MaterialProperty<Real> & _drho_g_dp;
  // The first derivative of liquid density wrt pressure
  MaterialProperty<Real> & _drho_l_dp;
  // The second derivative of gas density wrt pressure
  MaterialProperty<Real> & _drho_g_dp_2;
  // The second derivative of liquid density wrt pressure
  MaterialProperty<Real> & _drho_l_dp_2;
  // void_fraction
  MaterialProperty<Real> & _alpha;
};

#endif /* MOSKITOFLUIDWELL2P_H */
