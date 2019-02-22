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

#include "MoskitoMassFlowRateCoupled.h"

registerMooseObject("MoskitoApp", MoskitoMassFlowRateCoupled);

template <>
InputParameters
validParams<MoskitoMassFlowRateCoupled>()
{
  InputParameters params = validParams<NodalBC>();
  params.addRequiredParam<Real>("mass_flowrate",
        "The mass flowrate of the mixture (kg/s)");
  params.addRequiredParam<UserObjectName>("eos_uo",
        "The name of the userobject for EOS");
  params.addRequiredCoupledVar("pressure", "Pressure nonlinear variable (Pa)");
  params.addRequiredCoupledVar("enthalpy", "Specific enthalpy nonlinear variable (J/kg)");
  params.addClassDescription("Implements a NodalBC (Dirichlet) which calculates"
                            " mass weighted flow rate for momentum eq");
  return params;
}

MoskitoMassFlowRateCoupled::MoskitoMassFlowRateCoupled(const InputParameters & parameters)
  : NodalBC(parameters),
    eos_uo(getUserObject<MoskitoEOS2P>("eos_uo")),
    _m_dot(getParam<Real>("mass_flowrate")),
    _h(coupledValue("enthalpy")),
    _p(coupledValue("pressure")),
    _h_var_number(coupled("enthalpy")),
    _p_var_number(coupled("pressure"))
{
}

Real
MoskitoMassFlowRateCoupled::computeQpResidual()
{
  // Calculated mixture density
  Real One_by_rho_m;
  // Local variable of mass fraction
  Real vmfrac;
  // Temperature depending on p,H conditions
  Real T;

  eos_uo.VMFrac_from_p_h(_p[_qp], _h[_qp], vmfrac, T);
  One_by_rho_m = (vmfrac / eos_uo.rho_g_from_p_T(_p[_qp], T) + (1.0 - vmfrac) / eos_uo.rho_l_from_p_T(_p[_qp], T));

  return _u[_qp] - _m_dot * One_by_rho_m;
}

Real
MoskitoMassFlowRateCoupled::computeQpOffDiagJacobian(unsigned int jvar)
{
  // Local variable of mass fraction
  Real vmfrac;
  // Temperature depending on p,H conditions
  Real T;

  eos_uo.VMFrac_from_p_h(_p[_qp], _h[_qp], vmfrac, T);

  // The specific heat at constant pressure
  Real cp;
  cp = eos_uo.cp_m_from_p_T(_p[_qp], T, vmfrac);


  // Derivates of Densities
  Real rho_l, rho_g, drho_l_dp, drho_l_dT, drho_g_dp, drho_g_dT;

  eos_uo.rho_l_from_p_T(_p[_qp], T, rho_l, drho_l_dp, drho_l_dT);
  eos_uo.rho_g_from_p_T(_p[_qp], T, rho_g, drho_g_dp, drho_g_dT);

  // Assumption _mfrac has no derivative of p and h
  Real j = 0.0;
  if (jvar == _p_var_number)
    {
      j =  - vmfrac * drho_g_dp / (rho_g * rho_g);
      j -= (1.0 - vmfrac) * drho_l_dp / (rho_l * rho_l);
      j *= - _m_dot;
    }

  if (jvar == _h_var_number)
    {
      j =  - vmfrac * drho_g_dT / (rho_g * rho_g);
      j -= (1.0 - vmfrac) * drho_l_dT / (rho_l * rho_l);
      j *= - _m_dot / cp;
    }
    return j;
}
