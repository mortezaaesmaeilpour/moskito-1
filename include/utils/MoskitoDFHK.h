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

#ifndef MOSKITODFHK_H
#define MOSKITODFHK_H

#include "MoskitoDriftFlux.h"

class MoskitoDFHK : public MoskitoDriftFlux
{
public:
  MoskitoDFHK(const Real & v_m, const Real & rho_g, const Real & rho_l,
    const Real & mfrac, const Real & surf_ten, const Real & dia,
    const Real & dir, const Real & friction, const RealVectorValue & gravity,
    const RealVectorValue & well_dir);

  virtual void DFMCalculator() override;
  // Preparation of HK individual parameters
  void HKinitialisation(Real & v_sg, Real & v_sl, Real & vd_tb, Real & vd_b, Real & v_ms, Real & v_gb, Real & v_gc, Real & vd_mix);
  // final calculator
  void HKcalculator(const Real & vd_tb, const Real & v_sg, const Real & v_sl, const Real & v_ms, const Real & v_gc, const Real & v_gb, const Real & vd_b, const Real & vd_mix);
  // Calculation of Void fraction
  void HKvfrac(const Real & v_sg);

protected:
  // Main code determination of flow pattern, C0, vd and void fraction
  void cal_v_s(Real & v_sg, Real & v_sl);

  // Equations for calcualation of bubbly / Taylor rise velocities
  Real cal_vd_b();
  Real cal_vd_tb();
  Real cal_vd_mix(const Real & vd_b, const Real & vd_tb, const Real & v_gb, const Real & v_sg);

  // Calculation of thresholds of transitions
  Real cal_v_gb(const Real & v_sl, const Real & vd_b);
  Real cal_v_gc();
  Real cal_v_ms(const Real & v_sg);

  // Determine drift flow parameters for each flow pattern
  void Det_db_flow(const Real & vd_b);
  void Det_bubbly_flow(const Real & vd_b);
  void Det_churn_flow(const Real & v_ms, const Real & vd_mix);
  void Det_annular_flow(const Real & v_gc, const Real & v_sg);
  void Det_slug_flow(const Real & v_gb, const Real & v_sg, const Real & vd_mix);

  // Interpolation between different C0
  Real interpol(const Real & C0_1, const Real & C0_2, const Real & v_denom, const Real & v_num);

private:
  const Real _C0b = 1.2;
  const Real _C0db = 1.2;
  const Real _C0s_u = 1.2;
  const Real _C0s_d = 1.12;
  const Real _C0c_u = 1.15;
  const Real _C0c_d = 1.12;
  const Real _C0a = 1.0;

  Real _grav;
  // Angle between gravity vector and well_unity_vector
  Real _angle = 0.0; // fucntion to calculaute angle from gravity and well_unit_vect;
};

#endif
