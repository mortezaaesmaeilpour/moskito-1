/**************************************************************************/
/*  MOSKITO - Multiphysics cOupled Simulator toolKIT for wellbOres        */
/*                                                                        */
/*  Copyright (C) 2019 by Maziar Gholami Korzani                          */
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

#ifndef MOSKITOLATERALHEATXIONG_H
#define MOSKITOLATERALHEATXIONG_H

#include "Material.h"
#include "Function.h"

class MoskitoLateralHeatXiong;

template <>
InputParameters validParams<MoskitoLateralHeatXiong>();

class MoskitoLateralHeatXiong : public Material
{
public:
  MoskitoLateralHeatXiong(const InputParameters & parameters);
  virtual void computeQpProperties() override;
  Real Cal_ft(Real _alphaE, Real _rti, Real _lambdaE, Real _Uto, Real _rto);
  Real Cal_Te(Real Tsurf);

protected:
  // temperature
  const MaterialProperty<Real> & _T;
  // Radius tubing inner
  Real _rti;
  // Radius tubing outer
  MaterialProperty<Real> & _RadTubout;
  // Local parameter Radius tubing outer
  Real _rto;
  // Radius insulation
  Real _rins;
  // Radius casing inner
  Real _rci;
  // Radius casing outer respectivly cement
  Real _rcem;
  // Emissivity of inside casing surface
  Real _eao;
  // Emissivity of outside tubin/insulation surface
  Real _eai;
  // Density annulus fluid
  Real _rhoA;
  // Dyn Viscosity annulus fluid
  Real _nuA;
  // Thermal conductivity annulus fluid
  Real _lambdaA;
  // Heat capacity annulus fluid
  Real _cpA;
  //  Thermal volumetric expansion coefficient annulus fluid
  Real _betaA;
  // Annulus fluid thermal diffusivity = Temperaturleitfähigkeit
  Real _alphaA;
  // Thermal conductivoty formation (earth)
  Real _lambdaE;
  // Temperature at formation - cement boundary
  MaterialProperty<Real> & _Twb;
  // Formation thermal diffusivity
  Real _alphaE;
  // Thermal wellbore resistivity
  MaterialProperty<Real> & _Uto;
  //  Reziprok thermische Borhlochwiderstand
  MaterialProperty<Real> & _otU;
  // Surface temperature
  Real Tsurf;
  // Diameter of the pipe
  const MaterialProperty<Real> & _diameter;
  // Radius wellbore
  Real _rwb;
  // Thermal conductivity of the cement
  Real _lambdaCem;
  // Thermal conductivity of the casing
  Real _lambdaCas;
  // Thermal conductivity of the insulation
  Real _lambdaIns;
  // Thermal conductivity of the tubing
  Real _lambdaTub;
  // The gravity acceleration as a vector
  const MaterialProperty<RealVectorValue> & _gravity;
  // Function for geothermal gradient, has to be defined in the Input
  Function & gradT;
  // mixing approach
    MooseEnum _hc;
  enum HC_case {Dropkin_Sommerscales, Raithby_Hollands, Churchill};
// Ramey function parameter
  Real ft;

  // Convert from Si units to American system
  const Real m_to_ft    = 3.280839895;
  const Real m2_to_ft2  = 10.7639079;
  const Real m3_to_ft3  = 35.3146667;
  const Real kg_to_lbm  = 2.2046225;
  const Real gradC_to_gradR = 1.8;
  const Real Watt_to_Btu_per_h = 3.412141633;
  const Real J_to_Btu = 0.00094781712;
  const Real PI = 3.141592653589793238462643383279502884197169399375105820974944592308;
  const Real Rankine_absol = 491.67;
  // Stefan Bolzmann Konstante in SI units (W/(m² * K⁴))
  Real Boltz = 0.00000005670367;
  const Real s_to_h = 3600;
};

#endif
