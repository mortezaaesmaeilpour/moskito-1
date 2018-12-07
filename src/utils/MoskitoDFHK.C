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

MoskitoDFHK::MoskitoDFHK(const Real & v_m, const Real & rho_g,
  const Real & rho_l, const Real & mfrac, const Real & surf_ten, const Real & dia,
  const Real & dir, const Real & friction, const RealVectorValue & gravity,
  const RealVectorValue & well_dir)
  : MoskitoDriftFlux(v_m, rho_g, rho_l, mfrac, surf_ten, dia, dir, friction, gravity, well_dir),
    _grav(_gravity.norm())
{
}

void
MoskitoDFHK::DFMCalculator()
{
//
//   for (int i = 1; i < 6; ++i) //601
//     for (int j = 1; j < 5; ++j)
//     {
//       vm = 0.0003 * pow(10.0,(i / 1.0 ));    // kg/s = 0.198 x 10⁶ lbm/h
//       mf = 0.00001 * pow(10.0,(j / 1.0));
//
//
  // Superficial gas velocity in m/s
  Real v_sg= 0.0;
  // Superficial fluid velocity in m/s
  Real v_sl= 0.0;
  // velocity of dispersed bubbles
  Real vd_b= 0.0;
  // velocity of Taylor bubbles
  Real vd_tb= 0.0;
  // mixture velocity between Taylor bubbles velocity and dispersed bubble velocity for slug and churn flow
  Real vd_mix= 0.0;
  //  Thereshold for transition from bubbly to slug flow
  Real v_gb= 0.0;
  //  Thereshold for transition from churn to annular flow
  Real v_gc= 0.0;
  // Thereshold for transition from bubbly / slug flow to d_dubbly and slug to churn flow
  Real v_ms= 0.0;

  HKinitialisation(v_sg, v_sl, vd_tb, vd_b, v_ms, v_gb, v_gc, vd_mix);

  HKcalculator    (vd_tb, v_sg, v_sl, v_ms, v_gc, v_gb, vd_b, vd_mix);

  HKvfrac(v_sg);
//   std::cout<<vm * m_to_ft<<","<<mf<<","<<v_sg<<","<<v_sl<<","<<FlowPat<<","<<v_gc<<","<<v_ms<<","<<vd_b<<","<<vd_tb<<","<<v_gb<<std::endl;
// }
// abort();
}

void
MoskitoDFHK::HKinitialisation(Real & v_sg, Real & v_sl, Real & vd_tb, Real & vd_b, Real & v_ms, Real & v_gb, Real & v_gc, Real & vd_mix)
{
  cal_v_s(v_sg, v_sl);

  vd_tb = cal_vd_tb();
  vd_b = cal_vd_b();
  v_ms = cal_v_ms(v_sg);
  v_gc = cal_v_gc();
  v_gb = cal_v_gb(v_sl, vd_b);
  vd_mix = cal_vd_mix(vd_b, vd_tb, v_gb, v_sg);

  _vd = vd_b;
  _C0 = _C0b;

  HKvfrac(v_sg);
}

void
MoskitoDFHK::HKcalculator(const Real & vd_tb, const Real & v_sg, const Real & v_sl, const Real & v_ms, const Real & v_gc, const Real & v_gb, const Real & vd_b, const Real & vd_mix)
{
  _FlowPat = 0;

  // Check transition between slug and bubbly flow
  if (_vfrac <= 0.25 && vd_tb > _vd)
  {
      // Check transition between dispersed bubbly and bubbly flow
      if (v_ms < _v_m) //TODO - Check
        Det_db_flow(vd_b);
      else
        Det_bubbly_flow(vd_b);
  }
  else
  {
    // Check transition from slug to churn flow and from d_bubbly to churn flow
    if (v_sg > 1.08 * v_sl && v_ms < (_v_m * m_to_ft))
    {
      // Check transition between churn and annular flow
      if (v_sg < v_gc)
        Det_churn_flow(v_ms, vd_mix);
      else
        Det_annular_flow(v_gc, v_sg);
    }
    else
    {
      // Check transition between slug and dispersed bubbly flow
      if (_v_m < v_ms)
        Det_slug_flow(v_gb, v_sg, vd_mix);
      else
        Det_db_flow(vd_b);
    }
  }
}


void
MoskitoDFHK::HKvfrac(const Real & v_sg)
{
  _vfrac = v_sg / (_C0 * fabs(_v_m) + _vd * _dir);

  if (_vfrac < 0.0)
    _vfrac = 0.0;
}

void
MoskitoDFHK::cal_v_s(Real & v_sg, Real & v_sl)
{
  if (_dir == 1.0) // production, upflow
    v_sg =  _v_m * _rho_l * _mfrac  / ( - _mfrac * _rho_g + _rho_l * _mfrac + _rho_g);
  else                 // injection, downflow
    v_sg =  - _v_m * _rho_l * _mfrac / ( - _mfrac * _rho_g  - _rho_l * _mfrac + _rho_g);

  v_sg *= m_to_ft;

  if (_dir == 1.0) // production, upflow
    v_sl = _v_m * m_to_ft - v_sg;
  else          // injection, downflow
    v_sl = v_sg - _v_m * m_to_ft;
}

Real
MoskitoDFHK::cal_vd_b()
{
  return 1.53 * std::pow(m_to_ft * _grav * _surf_ten  * (_rho_l - _rho_g) / (std::pow(_rho_l,2.0) / m3_to_ft3),1.0 / 4.0);
}

Real
MoskitoDFHK::cal_vd_tb()
{
  return 0.35 * std::pow(_grav * m_to_ft * _dia * m_to_ft * (_rho_l - _rho_g) /
        _rho_l ,0.5) * std::pow(std::cos(_angle * M_PI / 180.0),0.5) *
        std::pow(1.0 + std::sin(_angle * M_PI / 180.0),1.2);
}

Real
MoskitoDFHK::cal_vd_mix(const Real & vd_b, const Real & vd_tb, const Real & v_gb, const Real & v_sg)
{
  return vd_b * (1.0 - std::exp(-0.1 * v_gb / (v_sg - v_gb))) +  vd_tb * std::exp(-0.1 * v_gb / (v_sg - v_gb));
}

Real
MoskitoDFHK::cal_v_gb(const Real & v_sl, const Real & vd_b)
{
  return ((_C0b / (4.0 - _C0b)) * v_sl + (1.0 / (4.0 - _C0b)) * vd_b * _dir) * cos(_angle * M_PI / 180.0);
}

Real
MoskitoDFHK::cal_v_gc()
{
  return 3.1 * std::pow(m_to_ft * _grav * _surf_ten * (_rho_l - _rho_g) / (std::pow(_rho_l,2.0) / m3_to_ft3),1.0 / 4.0);
}

Real
MoskitoDFHK::cal_v_ms(const Real & v_sg)
{
  return std::pow((0.725 + 4.15 * std::pow(v_sg / _v_m, 0.5)) /
          (2.0 * std::pow(0.4 * _surf_ten  / ((_rho_l - _rho_g ) *
          _grav / m3_to_ft3 * m_to_ft),0.5) * std::pow(_rho_l /
          (_surf_ten * m3_to_ft3) ,0.6) * std::pow(_friction /
          (2.0 * _dia * m_to_ft), 0.4)), 1.0 / 1.2);
}

void
MoskitoDFHK::Det_db_flow(const Real & vd_b)
{
  _FlowPat = 2;
  _C0 = _C0db;
  _vd = vd_b;
}

void
MoskitoDFHK::Det_bubbly_flow(const Real & vd_b)
{
  _FlowPat = 1;
  _C0 = _C0b;
  _vd = vd_b;
}

void
MoskitoDFHK::Det_churn_flow(const Real & v_ms, const Real & vd_mix)
{
  _FlowPat = 4;
  if (_dir == 1.0)
    _C0 = interpol(_C0s_u, _C0c_u, v_ms, _v_m);
  else
    _C0 = _C0c_d;

  _vd = vd_mix;
}

void
MoskitoDFHK::Det_annular_flow(const Real & v_gc, const Real & v_sg)
{
  _FlowPat = 5;
  if (_dir == 1.0)
    _C0 = interpol(_C0c_u, _C0a, v_gc, v_sg);
  else
    _C0 = interpol(_C0c_d, _C0a, v_gc, v_sg);

  _vd = 0.0;
}

void
MoskitoDFHK::Det_slug_flow(const Real & v_gb, const Real & v_sg, const Real & vd_mix)
{
  _FlowPat = 3;
  if (_dir == 1.0)
    _C0 = _C0s_u;
  else
    _C0 = interpol(_C0b, _C0s_d, v_gb, v_sg);

  _vd = vd_mix;
}

Real
MoskitoDFHK::interpol(const Real & C0_1, const Real & C0_2, const Real & v_denom, const Real & v_num)
{
  return C0_1 * (1 - std::exp(-0.1 * v_denom / (v_num - v_denom))) +  C0_2 * std::exp(-0.1 * v_denom / (v_num - v_denom));
}
