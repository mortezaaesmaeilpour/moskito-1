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

#ifndef MOSKITOMASS2P_H
#define MOSKITOMASS2P_H

#include "Kernel.h"

class MoskitoMass2P;

template <>
InputParameters validParams<MoskitoMass2P>();

class MoskitoMass2P : public Kernel
{
public:
  MoskitoMass2P(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned jvar) override;

  // The coupled flow_rate
  const VariableValue & _q_vol;
  // The coupled void_fraction
  const VariableValue & _alpha;

  // The gradient of the coupled flow_rate
  const VariableGradient & _grad_q_vol;
  // The gradient of the coupled void_fraction
  const VariableGradient & _grad_alpha;

  // Variable numberings
  unsigned _q_vol_var_number;
  unsigned _alpha_var_number;

  // The area of pipe
  const MaterialProperty<Real> & _area;
  // The unit vector of well direction
  const MaterialProperty<RealVectorValue> & _well_dir;
  // The density of gas
  const MaterialProperty<Real> & _rho_g;
  // The density of liquid
  const MaterialProperty<Real> & _rho_l;
  // The first derivative of gas density wrt pressure
  const MaterialProperty<Real> & _drho_g_dp;
  // The first derivative of liquid density wrt pressure
  const MaterialProperty<Real> & _drho_l_dp;
  // The second derivative of gas density wrt pressure
  const MaterialProperty<Real> & _drho_g_dp_2;
  // The second derivative of liquid density wrt pressure
  const MaterialProperty<Real> & _drho_l_dp_2;
};

#endif // MOSKITOMASS2P_H