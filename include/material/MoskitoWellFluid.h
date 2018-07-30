/**************************************************************************/
/*  TIGER - Hydro-thermal sImulator GEothermal Reservoirs                 */
/*                                                                        */
/*  Copyright (C) 2017 by Maziar Gholami Korzani                          */
/*  Karlsruhe Institute of Technology, Institute of Applied Geosciences   */
/*  Division of Geothermal Research                                       */
/*                                                                        */
/*  This file is part of TIGER App                                        */
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

#ifndef MOSKITOWELLFLUID_H
#define MOSKITOWELLFLUID_H

#define PI 3.141592653589793238462643383279502884197169399375105820974944592308

#include "Material.h"
#include "MoskitoEOS.h"

class MoskitoWellFluid;

template <>
InputParameters validParams<MoskitoWellFluid>();

class MoskitoWellFluid : public Material
{
public:
  MoskitoWellFluid(const InputParameters & parameters);
  virtual void computeQpProperties() override;

protected:
  /// Velocity in well
  MaterialProperty<Real> & _vel;
  /// Reynolds number in well
  MaterialProperty<Real> & _Re;
  /// Moody friction coefficient
  MaterialProperty<Real> & _friction;
  /// Well diameter
  MaterialProperty<Real> & _dia;
  /// Well area
  MaterialProperty<Real> & _area;

  /// Pressure and its fist and second derivatives
  MaterialProperty<Real> & _p;
  MaterialProperty<Real> & _dp_drho;
  MaterialProperty<Real> & _dp_dT;
  MaterialProperty<Real> & _dp_drho_2;
  MaterialProperty<Real> & _dp_dT_2;

  /// Userobject to equation of state
  const MoskitoEOS & _eos_UO;

  /// Density nonlinear variable value
  const VariableValue & _rho;
  /// Temperature nonlinear variable value
  const VariableValue & _T;
  /// Flow rate nonlinear variable value
  const VariableValue & _flow;

  /// function to calculate friction factor using Moody chart
  void MoodyFrictionFactor(Real & friction, Real rel_roughness, Real ReNo, MooseEnum roughness_type);

private:
  Real _d;
  Real _rel_roughness;
  MooseEnum _roughness_type;
};

#endif /* MOSKITOWELLFLUID_H */
