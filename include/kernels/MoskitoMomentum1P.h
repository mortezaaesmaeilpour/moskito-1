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

#ifndef MOSKITOMOMENTUM1P_H
#define MOSKITOMOMENTUM1P_H

#include "Kernel.h"

class MoskitoMomentum1P;

template <>
InputParameters validParams<MoskitoMomentum1P>();

class MoskitoMomentum1P : public Kernel
{
public:
  MoskitoMomentum1P(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned jvar) override;

  // The gradient of the coupled pressure
  const VariableGradient & _grad_p;
  // The gradient of the coupled specific enthalpy
  const VariableGradient & _grad_h;

  // Variable numberings
  unsigned _p_var_number;
  unsigned _h_var_number;

  // The specific heat at constant pressure
  const MaterialProperty<Real> & _cp;
  // The density
  const MaterialProperty<Real> & _rho;
  // The first derivative of density wrt pressure
  const MaterialProperty<Real> & _drho_dp;
  // The second derivative of density wrt pressure
  const MaterialProperty<Real> & _drho_dp_2;
  // The first derivative of density wrt temperature
  const MaterialProperty<Real> & _drho_dT;
  // The second derivative of density wrt temperature
  const MaterialProperty<Real> & _drho_dT_2;
  // The second derivative of density wrt temperature and pressure respectively
  const MaterialProperty<Real> & _drho_dTdp;
  // The pipe diameter
  const MaterialProperty<Real> & _d;
  // The pipe Moody friction factor
  const MaterialProperty<Real> & _f;
  // The gravity acceleration as a vector
  const MaterialProperty<RealVectorValue> & _gravity;
  // The area of pipe
  const MaterialProperty<Real> & _area;
  // The unit vector of well direction
  const MaterialProperty<RealVectorValue> & _well_dir;
};

#endif // MOSKITOMOMENTUM1P_H
