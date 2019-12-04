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

#include "MoskitoTemp2Enthalpy2P.h"

registerMooseObject("MoskitoApp", MoskitoTemp2Enthalpy2P);

template <>
InputParameters
validParams<MoskitoTemp2Enthalpy2P>()
{
  InputParameters params = validParams<InitialCondition>();
  // params.addRequiredParam<UserObjectName>("eos_uo",
  //       "The name of the userobject for EOS");
  // params.addRequiredParam<Real>("temperature", "Temperature gradient value of the IC");
  params.addRequiredParam<FunctionName>("geothermal_gradient",
        "Temperature gradient value of the IC");
  // params.addRequiredParam<Real>("mass_fraction", "Gradient of the Gas to Liquid ratio; value between 0 and 1"
                                // "0.0 = Liquid phase + critical; 1.0 = vapour phase ");
  params.addRequiredParam<FunctionName>("vmfrac_gradient", "Gradient of the Gas to Liquid ratio; value between 0 and 1"
                                "0.0 = Liquid phase + critical; 1.0 = vapour phase ");
  params.addRequiredCoupledVar("pressure", "Pressure nonlinear variable (Pa)");
  // params.declareControllable("temperature");
  params.addClassDescription("Implements a Initial Condition which calculates "
                            "specific enthalpy using temperature based on EOS "
                            " for 2 phase flow in Pure water");
  return params;
}

MoskitoTemp2Enthalpy2P::MoskitoTemp2Enthalpy2P(const InputParameters & parameters)
  : InitialCondition(parameters),
    _gradT(getFunction("geothermal_gradient")),
    _gradvmfrac(getFunction("vmfrac_gradient")),
    _p(coupledValue("pressure"))
{
  std::string eos_name = UserObjectName(name() + ":LiquidGas");
  {
    std::string class_name = "MoskitoWater97FluidProperties";
    InputParameters params = _app.getFactory().getValidParams(class_name);
    _fe_problem.addUserObject(class_name, eos_name, params);
  }
  _eos_uo = &_fe_problem.getUserObjectTempl<MoskitoWater97FluidProperties>(eos_name);
}

Real
MoskitoTemp2Enthalpy2P::value(const Point & _point)
{
  Real enthalpy = 0.0;
  if (_gradvmfrac.value(_t, _point) != 1.0 &&  _gradvmfrac.value(_t, _point) != 0.0)
  {
    std::cout<<"Temp1 = "<<_gradT.value(_t, _point)<<std::endl;
    Real pressure = _eos_uo->vaporPressure(_gradT.value(_t, _point));
    std::cout<<"vaporPressure = "<<pressure<<std::endl;
    Real hsatg = _eos_uo->h_from_p_T(pressure, _gradT.value(_t, _point), 2);
    Real hsatl = _eos_uo->h_from_p_T(pressure, _gradT.value(_t, _point), 1);
    enthalpy = hsatl + (hsatg - hsatl) * _gradvmfrac.value(_t, _point);
  }

  else
  {
    std::cout<<"Temp2 = "<<_gradT.value(_t, _point)<<std::endl;
    std::cout<<"presure = "<<_p[_qp]<<std::endl;
    unsigned int region = _eos_uo->inRegion(_p[_qp], _gradT.value(_t, _point));
    std::cout<<"region = "<<region<<std::endl;
    enthalpy = _eos_uo->h_from_p_T(_p[_qp], _gradT.value(_t, _point), region);
  }
  return enthalpy;
}
