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

#ifndef MOSKITOMOMENTUM2P_H
#define MOSKITOMOMENTUM2P_H

#include "Kernel.h"

class MoskitoMomentum2P;

template <>
InputParameters validParams<MoskitoMomentum2P>();

class MoskitoMomentum2P : public Kernel
{
public:
  MoskitoMomentum2P(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned jvar) override;

  // The coupled void_fraction
  const VariableValue & _alpha;

  // The gradient of the coupled pressure
  const VariableGradient & _grad_p;
  // The gradient of the coupled void_fraction
  const VariableGradient & _grad_alpha;

  // Variable numberings
  unsigned _p_var_number;
  unsigned _alpha_var_number;

  // The density of gas
  const MaterialProperty<Real> & _rho_g;
  // The density of liquid
  const MaterialProperty<Real> & _rho_l;
  // The first derivative of gas density wrt pressure
  const MaterialProperty<Real> & _drho_g_dp;
  // The first derivative of liquid density wrt pressure
  const MaterialProperty<Real> & _drho_l_dp;
  // The pipe diameter
  const MaterialProperty<Real> & _d;
  // The pipe Moody friction factor
  const MaterialProperty<Real> & _f;
  // The gravity acceleration as a vector
  RealVectorValue _gravity;
  // The area of pipe
  const MaterialProperty<Real> & _area;
  // The unit vector of well direction
  const MaterialProperty<RealVectorValue> & _well_dir;
};

#endif // MOSKITOMOMENTUM_H