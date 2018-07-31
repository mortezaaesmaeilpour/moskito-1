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

#ifndef MOSKITODENSITYCOUPLEDBC_H
#define MOSKITODENSITYCOUPLEDBC_H

#include "NodalBC.h"
#include "MoskitoEOS.h"
#include "Function.h"

class MoskitoDensityCoupledBC;

template <>
InputParameters validParams<MoskitoDensityCoupledBC>();

class MoskitoDensityCoupledBC : public NodalBC
{
public:
  MoskitoDensityCoupledBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  // Coupeled temperature
  const VariableValue & _T;

  // Variable numberings
  unsigned _T_var_number;

  // Pressure function
  Function & _p_func;

  // Userobject to equation of state
  const MoskitoEOS & _eos_UO;
};

#endif // MOSKITODENSITYCOUPLEDBC_H
