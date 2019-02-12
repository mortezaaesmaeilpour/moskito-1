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

#ifndef MOSKITOEOSNATURALGAS_H
#define MOSKITOEOSNATURALGAS_H

#include "MoskitoEOS1P.h"

class MoskitoEOSNaturalGas;

template <>
InputParameters validParams<MoskitoEOSNaturalGas>();

class MoskitoEOSNaturalGas : public MoskitoEOS1P
{
public:
  MoskitoEOSNaturalGas(const InputParameters & parameters);

  virtual Real rho_from_p_T(const Real & pressure, const Real & temperature) const override;
  virtual void rho_from_p_T(const Real & pressure, const Real & temperature,
                        Real & rho, Real & drho_dp, Real & drho_dT) const override;
  virtual Real T_to_h(const Real & temperature) const override;
  virtual Real h_to_T(const Real & enthalpy) const override;
  virtual Real cp(const Real & temperature) const override;
  virtual Real lambda(const Real & pressure, const Real & temperature) const override;
  void Pseudo_Critical_Calc(const Real & g);
  Real z_factor(const Real & pressure, const Real & temperature) const;


protected:
  // Molar mass of gas
  const Real _molar_mass;
  // Specific gravity
  const Real _gamma_g;
  // Universal Gas constant (J/mol.K)
  const Real _R;
  // Pseudo critical properties
  Real _P_pc;
  Real _T_pc;
  const Real _cp;
  const Real _lambda;

  // constants for z factor calculation based on Kareem et al 2016
  const std::array<Real, 20> a{
    {0.0 , 0.317842, 0.382216, -7.768354, 14.290531, 0.000002, -0.004693,
    0.096254, 0.166720, 0.966910, 0.063069, -1.966847, 21.0581, -27.0246,
    16.23, 207.783, -488.161, 176.29, 1.88453, 3.05921 }};
};

#endif /* MOSKITOEOSNATURALGAS_H */
