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

#ifndef MOSKITOCMOMENTUM_H
#define MOSKITOCMOMENTUM_H

#include "Kernel.h"

class MoskitoCMomentum;

template <>
InputParameters validParams<MoskitoCMomentum>();

class MoskitoCMomentum : public Kernel
{
public:
  MoskitoCMomentum(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned jvar) override;

  /// Coupled density (kg/m^3)
  const VariableValue & _rho;

  // Gradients of coupled temperature
  const VariableGradient & _grad_T;
  // Gradients of coupled density
  const VariableGradient & _grad_rho;

  // Variable numberings
  unsigned _T_var_number;
  unsigned _rho_var_number;

  // Well diameter
  const MaterialProperty<Real> & _d;
  // Well area
  const MaterialProperty<Real> & _area;
  // The first derivative of pressure wrt temperature
  const MaterialProperty<Real> & _dp_dT;
  // The first derivative of pressure wrt density
  const MaterialProperty<Real> & _dp_drho;
  // The second derivative of pressure wrt temperature
  const MaterialProperty<Real> & _dp_dT_2;
  // The second derivative of pressure wrt density
  const MaterialProperty<Real> & _dp_drho_2;
};

#endif // MOSKITOCMOMENTUM_H
