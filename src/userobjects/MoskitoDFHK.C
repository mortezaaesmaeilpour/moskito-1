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

#include "MoskitoDFHK.h"

registerMooseObject("MoskitoApp", MoskitoDFHK);

template <>
InputParameters
validParams<MoskitoDFHK>()
{
  InputParameters params = validParams<MoskitoDriftFlux>();
  params.addParam<Real>("surface_tension", 0.0,
        "Surface tension of the liquid phase (kg/s)");
  return params;
}

MoskitoDFHK::MoskitoDFHK(const InputParameters & parameters)
  : MoskitoDriftFlux(parameters),
    _surf_ten(getParam<Real>("surface_tension"))
{
  _surf_ten *= kg_to_lbm;
}

void
MoskitoDFHK::DFMCalculator(MoskitoDFGVar & input) const
{
  // conversion
  input._rho_g *= kg_to_lbm / m3m3_to_ft3;

  MoskitoHKLVar tmp;

  HKinitialisation(input, tmp);
  HKcalculator(input, tmp);
  HKvfrac(input, tmp);
  std::cout<<tmp.v_sg<<","<<tmp.v_sl<<","<<input._FlowPat<<","<<input._vfrac<<std::endl;
}

void
MoskitoDFHK::HKinitialisation(MoskitoDFGVar & input, MoskitoHKLVar & LVar) const
{
  cal_v_s(input, LVar);

  cal_vd_tb(input, LVar);
  cal_vd_b(input, LVar);
  cal_v_ms(input, LVar);
  cal_v_gc(input, LVar);
  cal_v_gb(input, LVar);
  cal_vd_mix(input, LVar);

  input._vd = LVar.vd_b;
  input._C0 = _C0b;

  HKvfrac(input, LVar);
}

void
MoskitoDFHK::HKcalculator(MoskitoDFGVar & input, MoskitoHKLVar & LVar) const
{
  input._FlowPat = 0;

  // Check transition between slug and bubbly flow
  if (input._vfrac <= 0.25 && LVar.vd_tb > input._vd)
  {
      // Check transition between dispersed bubbly and bubbly flow
      if (LVar.v_ms < (input._v_m * m_to_ft)) //TODO - Check
        Det_db_flow(input, LVar);
      else
        Det_bubbly_flow(input, LVar);
  }
  else
  {
    // Check transition from slug to churn flow and from d_bubbly to churn flow
    if (LVar.v_sg > 1.08 * LVar.v_sl && LVar.v_ms < (input._v_m * m_to_ft))
    {
      // Check transition between churn and annular flow
      if (LVar.v_sg < LVar.v_gc)
        Det_churn_flow(input, LVar);
      else
        Det_annular_flow(input, LVar);
    }
    else
    {
      // Check transition between slug and dispersed bubbly flow
      if ((input._v_m * m_to_ft) < LVar.v_ms)
        Det_slug_flow(input, LVar);
      else
        Det_db_flow(input, LVar);
    }
  }
}

void
MoskitoDFHK::HKvfrac(MoskitoDFGVar & input, MoskitoHKLVar & LVar) const
{
  input._vfrac = LVar.v_sg / (input._C0 * fabs(input._v_m * m_to_ft) + input._vd * input._dir);

  if (input._vfrac < 0.0)
    input._vfrac = 0.0;
}

void
MoskitoDFHK::cal_v_s(MoskitoDFGVar & input, MoskitoHKLVar & LVar) const
{
  if (input._dir == 1.0) // production, upflow
    LVar.v_sg =  input._v_m * input._rho_l * input._mfrac  / ( - input._mfrac * input._rho_g + input._rho_l * input._mfrac + input._rho_g);
  else                 // injection, downflow
    LVar.v_sg =  - input._v_m * input._rho_l * input._mfrac / ( - input._mfrac * input._rho_g  - input._rho_l * input._mfrac + input._rho_g);

  LVar.v_sg *= m_to_ft;

  if (input._dir == 1.0) // production, upflow
    LVar.v_sl = input._v_m * m_to_ft - LVar.v_sg;
  else          // injection, downflow
    LVar.v_sl = LVar.v_sg - input._v_m * m_to_ft;
}

void
MoskitoDFHK::cal_vd_b(MoskitoDFGVar & input, MoskitoHKLVar & LVar) const
{
  LVar.vd_b = 1.53 * std::pow(m_to_ft * input._grav * _surf_ten  * (input._rho_l - input._rho_g) / (std::pow(input._rho_l,2.0) / m3_to_ft3),1.0 / 4.0);
}

void
MoskitoDFHK::cal_vd_tb(MoskitoDFGVar & input, MoskitoHKLVar & LVar) const
{
  LVar.vd_tb = 0.35 * std::pow(input._grav * m_to_ft * input._dia * m_to_ft * (input._rho_l - input._rho_g) /
        input._rho_l ,0.5) * std::pow(std::cos(input._angle),0.5) *
        std::pow(1.0 + std::sin(input._angle),1.2);
}

void
MoskitoDFHK::cal_vd_mix(MoskitoDFGVar & input, MoskitoHKLVar & LVar) const
{
  LVar.vd_mix = LVar.vd_b * (1.0 - std::exp(-0.1 * LVar.v_gb / (LVar.v_sg - LVar.v_gb))) +  LVar.vd_tb * std::exp(-0.1 * LVar.v_gb / (LVar.v_sg - LVar.v_gb));
}

void
MoskitoDFHK::cal_v_gb(MoskitoDFGVar & input, MoskitoHKLVar & LVar) const
{
  LVar.v_gb = ((_C0b / (4.0 - _C0b)) * LVar.v_sl + (1.0 / (4.0 - _C0b)) * LVar.vd_b * input._dir) * cos(input._angle);
}

void
MoskitoDFHK::cal_v_gc(MoskitoDFGVar & input, MoskitoHKLVar & LVar) const
{
  LVar.v_gc = 3.1 * std::pow(m_to_ft * input._grav * _surf_ten * kg_to_lbm * (input._rho_l - input._rho_g) * kg_to_lbm / m3_to_ft3 / (std::pow(input._rho_g * kg_to_lbm / m3_to_ft3,2.0)),1.0 / 4.0);
}

void
MoskitoDFHK::cal_v_ms(MoskitoDFGVar & input, MoskitoHKLVar & LVar) const
{
  LVar.v_ms = std::pow((0.725 + 4.15 * std::pow(LVar.v_sg / input._v_m / m_to_ft, 0.5)) /
          (2.0 * std::pow(0.4 * _surf_ten  / ((input._rho_l - input._rho_g ) *
          input._grav / m3_to_ft3 * m_to_ft),0.5) * std::pow(input._rho_l /
          (_surf_ten * m3_to_ft3) ,0.6) * std::pow(input._friction /
          (2.0 * input._dia * m_to_ft), 0.4)), 1.0 / 1.2);
}

void
MoskitoDFHK::Det_db_flow(MoskitoDFGVar & input, MoskitoHKLVar & LVar) const
{
  input._FlowPat = 2;
  input._C0 = _C0db;
  input._vd = LVar.vd_b;
}

void
MoskitoDFHK::Det_bubbly_flow(MoskitoDFGVar & input, MoskitoHKLVar & LVar) const
{
  input._FlowPat = 1;
  input._C0 = _C0b;
  input._vd = LVar.vd_b;
}

void
MoskitoDFHK::Det_churn_flow(MoskitoDFGVar & input, MoskitoHKLVar & LVar) const
{
  input._FlowPat = 4;
  if (input._dir == 1.0)
    input._C0 = interpol(_C0s_u, _C0c_u, LVar.v_ms, input._v_m);
  else
    input._C0 = _C0c_d;

  input._vd = LVar.vd_mix;
}

void
MoskitoDFHK::Det_annular_flow(MoskitoDFGVar & input, MoskitoHKLVar & LVar) const
{
  input._FlowPat = 5;
  if (input._dir == 1.0)
    input._C0 = interpol(_C0c_u, _C0a, LVar.v_gc, LVar.v_sg);
  else
    input._C0 = interpol(_C0c_d, _C0a, LVar.v_gc, LVar.v_sg);

  input._vd = 0.0;
}

void
MoskitoDFHK::Det_slug_flow(MoskitoDFGVar & input, MoskitoHKLVar & LVar) const
{
  input._FlowPat = 3;
  if (input._dir == 1.0)
    input._C0 = _C0s_u;
  else
    input._C0 = interpol(_C0b, _C0s_d, LVar.v_gb, LVar.v_sg);

  input._vd = LVar.vd_mix;
}

Real
MoskitoDFHK::interpol(const Real & C0_1, const Real & C0_2, const Real & v_denom, const Real & v_num) const
{
  return C0_1 * (1 - std::exp(-0.1 * v_denom / (v_num - v_denom))) +  C0_2 * std::exp(-0.1 * v_denom / (v_num - v_denom));
}
