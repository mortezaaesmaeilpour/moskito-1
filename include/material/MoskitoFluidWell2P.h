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
#include "MoskitoEOS2P.h"

class MoskitoFluidWell2P;

template <>
InputParameters validParams<MoskitoFluidWell2P>();

class MoskitoFluidWell2P : public MoskitoFluidWellGeneral
{
public:
  MoskitoFluidWell2P(const InputParameters & parameters);
  virtual void computeQpProperties() override;

protected:
  // Userobject to equation of state
  const MoskitoEOS2P & eos_uo;
  // Userobject to Viscosity Eq
  // const MoskitoViscosity & viscosity_uo;

  // The specific heat of mixture at constant pressure
  MaterialProperty<Real> & _cp_m;
  // density of mixture
  MaterialProperty<Real> & _rho_m;
  // The first derivative of mixture density wrt pressure
  MaterialProperty<Real> & _drho_m_dp;
  // The second derivative of mixture density wrt pressure
  MaterialProperty<Real> & _drho_m_dp_2;
  // The first derivative of mixture density wrt temperature
  MaterialProperty<Real> & _drho_m_dT;
  // The second derivative of mixture density wrt temperature
  MaterialProperty<Real> & _drho_m_dT_2;
  // The second derivative of mixture density wrt temperature & pressure
  MaterialProperty<Real> & _drho_m_dTdp;
  // The second derivative of mixture density wrt pressure & temperature
  MaterialProperty<Real> & _drho_m_dpdT;
  // void_fraction
  MaterialProperty<Real> & _alpha;
private:
  // Gas density related properties
  Real _rho_g;
  Real _drho_g_dp;
  Real _drho_g_dp_2;
  Real _drho_g_dT;
  Real _drho_g_dT_2;
  Real _drho_g_dTdp;
  Real _drho_g_dpdT;
  // Liquid density related properties
  Real _rho_l;
  Real _drho_l_dp;
  Real _drho_l_dp_2;
  Real _drho_l_dT;
  Real _drho_l_dT_2;
  Real _drho_l_dTdp;
  Real _drho_l_dpdT;
};

#endif /* MOSKITOFLUIDWELL2P_H */
