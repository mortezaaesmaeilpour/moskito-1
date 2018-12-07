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

#ifndef MOSKITODRIFTFLUX_H
#define MOSKITODRIFTFLUX_H

#include "MoskitoFluidWell2P.h"

class MoskitoDriftFlux
{
public:
  MoskitoDriftFlux(const Real & v_m, const Real & rho_g, const Real & rho_l,
    const Real & mfrac, const Real & surf_ten, const Real & dia,
    const Real & dir, const Real & friction, const RealVectorValue & gravity,
    const RealVectorValue & well_dir);

  void DFMOutput(int & FlowPat, Real & vfrac, Real & C0, Real & vd);
  virtual void DFMCalculator() = 0;

  // Flow pattern 0 = nothing, 1 = bubbly, 2 = dispersed_bubbly, 3 = slug, 4 = churn, 5 = annular
  int _FlowPat;
  // Volumetric fraction of void phase
  Real _vfrac;
  // Drift Flux parameters
  Real _C0;
  Real _vd;

protected:
  // mixture velocity of 2P-system
  const Real _v_m;
  // Gas density
  const Real _rho_g;
  // Liquid density
  const Real _rho_l;
  // Mass fraction of void phase
  const Real _mfrac;
  // Surface tension of liquid phase
  const Real _surf_ten;
  // Well diameter
  const Real _dia;
  // Flow direction
  const Real _dir;
  // Well friction
  const Real _friction;
  // The gravity acceleration as a vector
  const RealVectorValue _gravity;
  // unit vector along well
  const RealVectorValue _well_dir;

  // Convert from Si units to American system
  const Real m_to_ft    = 3.2808398;
  const Real m3_to_ft3  = 35.3146667;
  const Real kg_to_lbm  = 2.2046225;
};

#endif
