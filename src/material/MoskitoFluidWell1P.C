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

#include "MoskitoFluidWell1P.h"

registerMooseObject("MoskitoApp", MoskitoFluidWell1P);

template <>
InputParameters
validParams<MoskitoFluidWell1P>()
{
  InputParameters params = validParams<Material>();

  params.addRequiredCoupledVar("pressure", "Pressure nonlinear variable (Pa)");
  params.addRequiredCoupledVar("flow_rate", "Mixture flow rate nonlinear variable (m^3/s)");
  params.addCoupledVar("temperature", 273.15, "Temperature nonlinear variable (K)");

  params.addRequiredRangeCheckedParam<Real>("well_diameter", "well_diameter>0", "Well diameter (m)");
  params.addRangeCheckedParam<Real>("roughness", 2.5e-5, "roughness>0", "Material roughness of well casing (m)");
  params.addRangeCheckedParam<Real>("manual_friction_factor", 0.0, "manual_friction_factor>=0",
                                    "User defined constant friction factor (if it is defined, the automatic "
                                    " moody friction factor based on roughness and type of casing will be disabled)");
  params.addRequiredParam<UserObjectName>("eos_UO", "The name of the userobject for EOS");
  params.addRequiredParam<UserObjectName>("viscosity_UO",
                                          "The name of the userobject for Viscosity Eq");

  MooseEnum RT("rough=1 smooth=2");
  params.addParam<MooseEnum>("roughness_type", RT="smooth", "Well casing roughness type [rough, smooth].");

  MooseEnum WD("x=1 -x=2 y=3 -y=4 z=5 -z=6");
  params.addRequiredParam<MooseEnum>(
      "well_direction", WD, "Well dominent direction towards bottom hole [x, -x, y, -y, z, -z].");

  return params;
}

MoskitoFluidWell1P::MoskitoFluidWell1P(const InputParameters & parameters)
  : Material(parameters),
    _vel(declareProperty<Real>("well_velocity")),
    _Re(declareProperty<Real>("well_reynolds_no")),
    _friction(declareProperty<Real>("well_moody_friction")),
    _dia(declareProperty<Real>("well_diameter")),
    _area(declareProperty<Real>("well_area")),
    _rho(declareProperty<Real>("density")),
    _drho_dp(declareProperty<Real>("drho_dp")),
    _drho_dp_2(declareProperty<Real>("drho_dp_2")),
    _well_unit_vect(declareProperty<RealVectorValue>("well_direction_vector")),
    _eos_UO(getUserObject<MoskitoEOS>("eos_UO")),
    _viscosity_UO(getUserObject<MoskitoViscosity>("viscosity_UO")),
    _T(coupledValue("temperature")),
    _P(coupledValue("pressure")),
    _flow(coupledValue("flow_rate")),
    _d(getParam<Real>("well_diameter")),
    _rel_roughness(getParam<Real>("roughness")),
    _f(getParam<Real>("manual_friction_factor")),
    _f_defined(parameters.isParamSetByUser("manual_friction_factor")),
    _roughness_type(getParam<MooseEnum>("roughness_type")),
    _well_direction(getParam<MooseEnum>("well_direction"))
{
  _rel_roughness /= _d;
}

void
MoskitoFluidWell1P::computeQpProperties()
{
  Real temp;

  _eos_UO.drho_dpT(_P[_qp], _T[_qp], _rho[_qp], _drho_dp[_qp], temp);
  _eos_UO.drho_dpT_2(_P[_qp], _T[_qp], _drho_dp_2[_qp], temp);

  _dia[_qp] = _d;
  _area[_qp] = PI * _d * _d / 4.0;

  _vel[_qp] = _flow[_qp] / _area[_qp];
  _Re[_qp] = _rho[_qp] * _dia[_qp] * fabs(_vel[_qp]) / _viscosity_UO.mu(_P[_qp], _T[_qp]);
  if (_f_defined)
    _friction[_qp] = _f;
  else
    MoodyFrictionFactor(_friction[_qp], _rel_roughness, _Re[_qp], _roughness_type);

  _well_unit_vect[_qp] = WellUnitVector();
}

void
MoskitoFluidWell1P::MoodyFrictionFactor(Real & friction, Real rel_roughness, Real ReNo, MooseEnum roughness_type)
{
  if (ReNo > 0.0)
  {
    if (ReNo < 3500.0)
      friction = 64.0 / ReNo;
    else
      switch (roughness_type)
        {
          case 1:
            Real a, b, c, d;
            a = -2.0 * std::log10(rel_roughness / 3.7 + 12.0 / ReNo);
            b = -2.0 * std::log10(rel_roughness / 3.7 + 2.51 * a / ReNo);
            c = -2.0 * std::log10(rel_roughness / 3.7 + 2.51 * b / ReNo);
            d = a - std::pow(b - a,2.0) / (c - 2.0 * b + a);
            friction = std::pow(1.0 / d,2.0);
            break;

          case 2:
            friction = 0.184 * std::pow(ReNo,-0.2);
            break;
        }
  }
  else
    friction = 0.0;
}

RealVectorValue
MoskitoFluidWell1P::WellUnitVector()
{
  RealVectorValue p0, p1, p;
  p0 = _current_elem->point(0);
  p1 = _current_elem->point(1);

  switch (_well_direction)
  {
    case 1:
      if (p0(0) > p1(0))
        p = p0 - p1;
      else
        p = p1 - p0;
      break;

    case 2:
      if (p0(0) < p1(0))
        p = p0 - p1;
      else
        p = p1 - p0;
      break;

    case 3:
      if (p0(1) > p1(1))
        p = p0 - p1;
      else
        p = p1 - p0;
      break;

    case 4:
      if (p0(1) < p1(1))
        p = p0 - p1;
      else
        p = p1 - p0;
      break;

    case 5:
      if (p0(2) > p1(2))
        p = p0 - p1;
      else
        p = p1 - p0;
      break;

    case 6:
      if (p0(2) < p1(2))
        p = p0 - p1;
      else
        p = p1 - p0;
      break;
  }

  p /= p.norm();
  return p;
}