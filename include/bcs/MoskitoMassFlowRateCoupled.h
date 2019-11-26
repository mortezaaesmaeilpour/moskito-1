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

#ifndef MOSKITOMASSFLOWRATECOUPLED_H
#define MOSKITOMASSFLOWRATECOUPLED_H

#include "NodalBC.h"
#include "MoskitoEOS2P.h"

class MoskitoMassFlowRateCoupled;

template <>
InputParameters validParams<MoskitoMassFlowRateCoupled>();

class MoskitoMassFlowRateCoupled : public NodalBC
{
public:
  MoskitoMassFlowRateCoupled(const InputParameters & parameters);
  virtual Real computeQpResidual() override;
  // virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

protected:
  // Userobject to equation of state
  const MoskitoEOS2P & eos_uo;
  //Mass flow rate assigned by the user
  const Real & _m_dot;
  // Reading of coupled parameters
  const VariableValue & _h;
  const VariableValue & _p;

};

#endif // MOSKITOMASSFLOWRATE_H
