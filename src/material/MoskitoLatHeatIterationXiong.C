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

#include "MoskitoLatHeatIterationXiong.h"
#include "Conversion.h"

registerMooseObject("MoskitoApp", MoskitoLatHeatIterationXiong);

template <>
InputParameters
validParams<MoskitoLatHeatIterationXiong>()
{
    InputParameters params = validParams<Material>();
    params += validParams<NewtonIteration>();
    params.addClassDescription("Materials for the Lateral heat transfer between "
          "wellbore and formation");
    params.addRequiredParam<Real>("radius_tubbing_inner",
          "Inner radius of the tubing (m)");
    params.addRequiredParam<Real>("radius_tubbing_outer",
          "Outer radius of the tubing (m)");
    params.addParam<Real>("radius_insulation", 0.0,
          "Radius of the insulation (m), if not existing put 0.0");
    params.addParam<Real>("radius_casing_inner", 0.0,
          "Inner radius of the casing (m)");
    params.addParam<Real>("radius_cement", 0.0,
          "Outer radius of the casing respectivly inner radius cement (m)");
    params.addParam<Real>("emissivity_annulus_outer", 0.9,
          "Emissivity of inside casing surface ()");
    params.addParam<Real>("emissivity_annulus_inner", 0.9,
          "Emissivity of outside tubing/insulation surface ()");
    params.addParam<Real>("density_annulus", 1000,
          "Density of the formation (kg/m³))");
    params.addParam<Real>("dyn_viscosity_annulus",0.0002843,
          "Dynamische viscosity of the annulus fluid (kg/(m*s))");
    params.addParam<Real>("conductivity_annulus", 0.6 ,
          "Thermal conductivity of the annulus fluid (W/(m*K))");
    params.addParam<Real>("capacity_annulus", 4181.3,
          "Specific Heat capacity of the annulus fluid (J/(K*kg))");
    params.addParam<Real>("thermal_expansion_annulus", 0.000207,
          "Thermal volumetric expansion coefficient annulus fluid (1 / K )");
    params.addParam<Real>("thermal_diffusivity_annulus", 0.00000014,
          "Thermal diffusivity of the annulus fluid (m²/s)");
    params.addRequiredParam<Real>("conductivity_rock",
          "Thermal conductivity of the formation (W/(m*K))");
    params.addRequiredParam<Real>("Surface_temperature",
          "Surface temperature (°C); compare with BC");
    params.addRequiredParam<Real>("thermal_diffusivity_rock",
          "Thermal diffusivity of the formation (m²/s)");
    params.addParam<Real>("conductivity_cement",46.5,
          "Thermal conductivity of the cement (W/(m*K))");
    params.addParam<Real>("conductivity_casing",46.5,
          "Thermal conductivity of the casing (W/(m*K))");
    params.addParam<Real>("conductivity_insulation",0.5,
          "Thermal conductivity of the insulation (W/(m*K))");
    params.addRequiredParam<Real>("conductivity_tubing",
          "Thermal conductivity of the tubing (W/(m*K))");
    params.addRequiredParam<FunctionName>("geothermal_gradient",
          "a function to define the geothermal gradient (°C/100m)");
    MooseEnum hc_model
          ("Dropkin_Sommerscales Raithby_Hollands Churchill");
    params.addParam<MooseEnum>("hc_calucation_model", hc_model,
          "Type of calculation of hc; Methodology: Dropkin & Sommerscales 1965,"
          "Raithby & Holland 1975, Churchill 1983; While DS and RH are similar in the results, Churchill approach differs strongly");
    MooseEnum time_type
          ("user_time simulation_time");
    params.addParam<MooseEnum>("time_model", time_type,
          "Define time for steady simulations"
          "user_time, simulation_time");
    params.addParam<Real>("user_defined_time", 86400,
          "Time defined by the user for steady state simulation, Default = 1 day");
    params.addParam<Real>("derivative_tolerance", 0.001,
          "Tolerance to calculate derivatives based on numerical differentiation");
    return params;
}

MoskitoLatHeatIterationXiong::MoskitoLatHeatIterationXiong(const InputParameters & parameters)
  : Material(parameters),
    NewtonIteration(parameters),
    _T(getMaterialProperty<Real>("temperature")),
    _rti(getParam<Real>("radius_tubbing_inner")),
    _RadTubout(declareProperty<Real>("radius_tubbing_outer")),
    _TRankine(declareProperty<Real>("Temperature_Rankine")),
    _rto(getParam<Real>("radius_tubbing_outer")),
    _rins(getParam<Real>("radius_insulation")),
    _rci(getParam<Real>("radius_casing_inner")),
    _rcem(getParam<Real>("radius_cement")),
    _eao(getParam<Real>("emissivity_annulus_outer")),
    _eai(getParam<Real>("emissivity_annulus_inner")),
    _rhoAnnulus(getParam<Real>("density_annulus")),
    _nuAnnulus(getParam<Real>("dyn_viscosity_annulus")),
    _lambdaAnnulus(getParam<Real>("conductivity_annulus")),
    _cpAnnulus(getParam<Real>("capacity_annulus")),
    _betaAnnulus(getParam<Real>("thermal_expansion_annulus")),
    _alphaAnnulus(getParam<Real>("thermal_diffusivity_annulus")),
    _lambdaRock(getParam<Real>("conductivity_rock")),
    _alphaRock(getParam<Real>("thermal_diffusivity_rock")),
    _Uto(declareProperty<Real>("thermal_resistivity_well")),
    _Twb(declareProperty<Real>("temperature_well_formation_interface")),
    _Tsurf(getParam<Real>("Surface_temperature")),
    _diameter(getMaterialProperty<Real>("well_diameter")),
    _lambdaCem(getParam<Real>("conductivity_cement")),
    _lambdaCas(getParam<Real>("conductivity_casing")),
    _lambdaIns(getParam<Real>("conductivity_insulation")),
    _lambdaTub(getParam<Real>("conductivity_tubing")),
    _gravity(getMaterialProperty<RealVectorValue>("gravity")),
    gradT(getFunction("geothermal_gradient")),
    _hc(getParam<MooseEnum>("hc_calucation_model")),
    _td(getParam<MooseEnum>("time_model")),
    _ut(getParam<Real>("user_defined_time")),
    _tol(getParam<Real>("derivative_tolerance"))
{
  _Tsurf *= gradC_to_gradR;
  _Tsurf += Rankine_absol;
  Boltz = Boltz * Watt_to_Btu_per_h / (m2_to_ft2 * gradC_to_gradR * gradC_to_gradR * gradC_to_gradR * gradC_to_gradR);
  _rti *= m_to_ft;
  _rto *= m_to_ft;
  _rins *= m_to_ft;
  _rci *= m_to_ft;
  _rcem *= m_to_ft;
  _lambdaRock *= Watt_to_Btu_per_h / (m_to_ft * gradC_to_gradR);
  _alphaRock *= m2_to_ft2;
  _rhoAnnulus *= kg_to_lbm / m3_to_ft3;
  _nuAnnulus *=  kg_to_lbm / m_to_ft * s_to_h;
  _lambdaAnnulus *=  Watt_to_Btu_per_h / (m_to_ft * gradC_to_gradR);
  _cpAnnulus *=  J_to_Btu / (kg_to_lbm * gradC_to_gradR);
  _betaAnnulus /= gradC_to_gradR;
  _alphaAnnulus *= m2_to_ft2 * s_to_h;
  _lambdaCem *=  Watt_to_Btu_per_h / (m_to_ft * gradC_to_gradR);
  _lambdaCas *=  Watt_to_Btu_per_h / (m_to_ft * gradC_to_gradR);
  _lambdaIns *=  Watt_to_Btu_per_h / (m_to_ft * gradC_to_gradR);
  _lambdaTub *=  Watt_to_Btu_per_h / (m_to_ft * gradC_to_gradR);
}

// Compute dimensionless time (ft)
Real
MoskitoLatHeatIterationXiong::transienttimefunction(Real Uto)
{
  Real Timing;
  if (_td == user_time)
    Timing = _ut;
  else
    Timing = _t;

  if (Timing < 604800)
  {
  Real column;
  if(Uto == 0.0) // Case can only occur in the initial step
    column = _rto / _lambdaRock;
  else
    column = _rto * Uto / _lambdaRock;

  Real row = _alphaRock * Timing / (_rti * _rti);

  int col_number = 0;
  int row_number = 0;

  const std::array<Real, 14> column_vec{
    {0.01, 0.02, 0.05,
     0.1,  0.2,  0.5,
     1.0,  2.0,  5.0,
     10.0, 20.0, 50.0,
     100.0, 1000.0}};

  const std::array<Real, 10> row_vec{
    {0.1,  0.2,  0.5,
     1.0,  2.0,  5.0,
     10.0, 20.0, 50.0,
     100.0}};

  const std::array<std::array<Real, 14>, 10> Ramey{
   {{{0.313, 0.313, 0.314, 0.316, 0.318, 0.323, 0.330, 0.345, 0.373, 0.396, 0.417, 0.433, 0.438, 0.445}},
    {{0.423, 0.423, 0.424, 0.427, 0.430, 0.439, 0.452, 0.473, 0.511, 0.538, 0.568, 0.572, 0.578, 0.588}},
    {{0.616, 0.617, 0.619, 0.623, 0.629, 0.644, 0.666, 0.698, 0.745, 0.772, 0.790, 0.802, 0.806, 0.811}},
    {{0.802, 0.803, 0.806, 0.811, 0.820, 0.842, 0.872, 0.910, 0.958, 0.984,  1.00,  1.01,  1.01,  1.02}},
    {{ 1.02,  1.02,  1.03,  1.04,  1.05,  1.08,  1.11,  1.15,  1.20,  1.22,  1.24,  1.24,  1.25,  1.25}},
    {{ 1.36,  1.37,  1.37,  1.38,  1.40,  1.44,  1.48,  1.52,  1.56,  1.57,  1.58,  1.59,  1.59,  1.59}},
    {{ 1.65,  1.66,  1.66,  1.67,  1.69,  1.73,  1.77,  1.81,  1.84,  1.86,  1.86,  1.87,  1.87,  1.88}},
    {{ 1.96,  1.97,  1.97,  1.99,  2.00,  2.05,  2.09,  2.12,  2.15,  2.16,  2.16,  2.17,  2.17,  2.17}},
    {{ 2.39,  2.39,  2.40,  2.42,  2.44,  2.48,  2.51,  2.54,  2.56,  2.57,  2.57,  2.57,  2.58,  2.58}},
    {{ 2.73,  2.73,  2.74,  2.75,  2.77,  2.81,  2.84,  2.86,  2.88,  2.89,  2.89,  2.89,  2.89,  2.90}}}};

  for(int i = 0; i <= 12; ++i){
    if (column > column_vec[13]){
      col_number = 13;
      break;
    }
    else if (column > column_vec[i])
      col_number = i;
    else
      break;
    }

  for(int j = 0; j <= 8; ++j){
    if (row > row_vec[9]){
      row_number = 9;
      break;
    }
    else if (row > row_vec[j])
      row_number = j;
    else
      break;
    }

   // Horizontal interpolation 1
   Real Hor1;
   if (col_number == 13)
     Hor1 = Ramey[row_number][col_number];
   else {
     Hor1 = (std::log(column) -  std::log(column_vec[col_number]));
     Hor1 /= (std::log(column_vec[col_number+1]) - std::log(column_vec[col_number]));
     Hor1 *= (Ramey[row_number][col_number+1] - Ramey[row_number][col_number]);
     Hor1 += Ramey[row_number][col_number];
   }

   // Horizontal interpolation 2
   Real Hor2;
   if (col_number == 13)
     Hor2 = Ramey[row_number+1][col_number];
   else {
     Hor2 = (std::log(column) -  std::log(column_vec[col_number]));
     Hor2 /= (std::log(column_vec[col_number+1]) - std::log(column_vec[col_number]));
     Hor2 *= (Ramey[row_number+1][col_number+1] - Ramey[row_number+1][col_number]);
     Hor2 += Ramey[row_number+1][col_number];
   }

   if (column < column_vec[0]){
    Hor1 = Ramey[row_number][col_number];
    Hor2 = Ramey[row_number+1][col_number];
   }

   if (row < row_vec[0])
     Hor2 = Hor1;

   if (row_number == 9)
     ft = Hor1;
   else {
     ft = (std::log(row) -  std::log(row_vec[row_number]));
     ft /= (std::log(row_vec[row_number+1]) - std::log(row_vec[row_number]));
     ft *= Hor2 - Hor1;
     ft += Hor1;
   }
  }
  else
    ft = std::log(2.0 * std::pow(_alphaRock * Timing,0.5) / _rwb) - 0.29;

  return ft;

}

Real
MoskitoLatHeatIterationXiong::TemperatureFormation(Real _Tsurf)
{
  Real Test;
  Test = _Tsurf +  gradT.value(_t, _q_point[_qp]) * gradC_to_gradR;
  return _Tsurf +  gradT.value(_t, _q_point[_qp]) * gradC_to_gradR;
}

Real
MoskitoLatHeatIterationXiong::TemperatureWBinterface(Real Uto, Real TRankine)
{
  Real Twb = 0.0;
  Twb += _rto * Uto * transienttimefunction(Uto) * TRankine;
  Twb += _lambdaRock * TemperatureFormation(_Tsurf);
  Twb /= _rto * Uto * transienttimefunction(Uto) + _lambdaRock;
  // std::cout<<"rto = "<< _rto<<std::endl;
  // std::cout<<"TRankine = "<< TRankine<<std::endl;
  // std::cout<<"_lambdaRock = "<< _lambdaRock<<std::endl;

  return Twb;
}

Real
MoskitoLatHeatIterationXiong::TemperatureCasingAnnulusInterface(Real Uto, Real TRankine)
{
  Real Tci = 0.0;
  if (_rcem != 0.0)
  {
    Tci += std::log(_rwb / _rcem) / _lambdaCem;
    Tci += std::log(_rcem / _rci) / _lambdaCas;
    Tci *= _rto * Uto * (TRankine - TemperatureWBinterface(Uto, TRankine));
  }
  Tci += TemperatureWBinterface(Uto, TRankine);
  return Tci;
}

Real
MoskitoLatHeatIterationXiong::TemperatureTubingOuter(Real Uto, Real TRankine)
{
  Real Tto;
  Tto = std::log(_rto / _rti) / _lambdaTub;
  Tto *= - _rto * Uto * (TRankine - TemperatureWBinterface(Uto, TRankine));
  Tto += TRankine;
  return Tto;
}

Real
MoskitoLatHeatIterationXiong::RadialHeatTransferCoefficient(Real Uto, Real TRankine)
{
  Real OverFtci = 0.0;
  OverFtci += 1.0 / _eao - 1.0;
  OverFtci *= _rai / _rao;
  OverFtci += 1.0 / _eai;

  Real hr;
  hr = Boltz / OverFtci;
  hr *= TemperatureTubingOuter(Uto, TRankine) + TemperatureCasingAnnulusInterface(Uto, TRankine);
  hr *= TemperatureTubingOuter(Uto, TRankine) * TemperatureTubingOuter(Uto, TRankine) + TemperatureCasingAnnulusInterface(Uto, TRankine) * TemperatureCasingAnnulusInterface(Uto, TRankine);
  return hr;
}

Real
MoskitoLatHeatIterationXiong::ConvectiveHeatTransferCoefficient(Real Uto, Real TRankine, Real grav)
{
  Real khc, hc;
  // Calculate Prandtl number
  Real Pr;
    Real Test1;
  Pr = _cpAnnulus * _nuAnnulus / _lambdaAnnulus;

  switch(_hc)
  {
  case HC_case::Dropkin_Sommerscales:

  khc = 0.049 * std::pow(Grashof(Uto, TRankine, grav) * Pr, 0.333);
  khc *= std::pow(Pr,0.074) * _lambdaAnnulus;
  hc = khc / (_rai * std::log(_rao / _rai));
  break;

  case HC_case::Raithby_Hollands:

  //  I am quite unclear about 0.384. In Xiong it is 0.386 but checking the original document it should be 0.48 * 0.8
  khc = 0.384 * std::pow(Pr / (0.861 + Pr),0.25);
  khc *= std::pow(Rayleigh(grav, Uto, TRankine),0.25) * _lambdaAnnulus;

  hc = khc / (_rai * std::log(_rao / _rai));
  break;

  // Document not available
  case HC_case::Churchill:
  Real fPR, Nu;
  fPR = std::pow(1 + std::pow(0.5 / Pr,9.0 / 16.0), - 16.0 / 9.0);
  // Calculate Nusselt number
  Nu = 0.364 * std::pow(Rayleigh(grav, Uto, TRankine) * fPR, 0.25);
  Nu *= std::pow(_rao / _rai,0.5);

  hc = Nu * _lambdaAnnulus / (2.0 * _rao);
  break;
  }
  return hc;
}

Real
MoskitoLatHeatIterationXiong::Grashof(Real Uto, Real TRankine, Real grav)
{
  Real Gr;
  Gr = (_rao - _rai) * (_rao - _rai) * (_rao - _rai);
  Gr *= _rhoAnnulus * _rhoAnnulus * _betaAnnulus;
  Gr *= TemperatureTubingOuter(Uto, TRankine) - TemperatureCasingAnnulusInterface(Uto, TRankine);
  Gr *= grav * m_to_ft * s_to_h * s_to_h;
  Gr /= _nuAnnulus * _nuAnnulus;
  Real TestA = _rao - _rai;
  return Gr;
}

Real
MoskitoLatHeatIterationXiong::Rayleigh(Real grav, Real Uto, Real TRankine)
{
  Real Lc, Ray;
  Lc = 2 * std::pow(std::log(_rao / _rai),4.0 / 3.0);
  Lc /= std::pow(std::pow(_rai,-3.0 / 5.0) + std::pow(_rao,-3.0 / 5.0), 5.0 / 3.0);

  Ray = grav * m_to_ft * s_to_h * s_to_h;
  Ray *= _betaAnnulus * (TemperatureTubingOuter(Uto, TRankine) - TemperatureCasingAnnulusInterface(Uto, TRankine)) * Lc * Lc * Lc;
  Ray /= _alphaAnnulus *  _nuAnnulus / _rhoAnnulus;
  return Ray;
}

Real
MoskitoLatHeatIterationXiong::computeReferenceResidual(const Real trail_value, const Real scalar)
{
  return 1e-2;
}

Real
MoskitoLatHeatIterationXiong::computeResidual(const Real trail_value, const Real scalar)
{
  Real hr;
  hr =  RadialHeatTransferCoefficient(scalar, _TRankine[_qp]);
  Real hc, grav;
  grav = _gravity[_qp].norm();
  hc = ConvectiveHeatTransferCoefficient(scalar, _TRankine[_qp], grav);

  // Auxillary variable
  Real Aux, Uto;
  // Calculation of thermal wellbore resistivity
  Aux =  std::log(_rto / _rti) / _lambdaTub;
  Aux += 1.0 / (_rai * (hc + hr));
  if (_rins != 0.0)
    Aux += std::log(_rins / _rto) / _lambdaIns;
  if (_rci != 0.0)
    Aux += std::log(_rcem / _rci) / _lambdaCas;
  if (_rcem != 0.0)
    Aux += std::log(_rwb / _rcem) / _lambdaCem;

  Uto = 1.0 / (Aux * _rto);

  // std::cout<<"hc = "<< hc<<std::endl;
  // std::cout<<"z coord ="<<_q_point[_qp] <<std::endl;
  // std::cout<<"Tto = "<< Tto<<std::endl;
  // std::cout<<"_T = "<< _T[_qp]<<std::endl;
  // std::cout<<"Uto = "<< Uto<<std::endl;
  // std::cout<<"TRankine = "<< TRankine<<std::endl;
  // std::cout<<"_lambdaRock = "<< _lambdaRock<<std::endl;
  // std::cout<<"hc = "<< hc<<std::endl;
  // std::cout<<"z coord ="<<_q_point[_qp] <<std::endl;
  // Real Test1 =  std::log(_rwb / _rcem) / _lambdaCem;;
  // Real Test2 = _rwb / _rcem;
  // Real Test3 = std::log(_rwb / _rcem);

  // std::cout<<"Tci = "<< Tci<<std::endl;
  // std::cout<<"Test1 = "<< Test1<<std::endl;
  // std::cout<<"LCem = "<< _lambdaCem<<std::endl;
  // std::cout<<"TRankine = "<< TRankine<<std::endl;
  // std::cout<<"Test2 = "<< Test2<<std::endl;
  // std::cout<<"Test3 = "<< Test3<<std::endl;
  // std::cout<<"_rcem = "<< _rcem<<std::endl;
  // std::cout<<"Twb = "<< Twb<<std::endl;
  return Uto - scalar;
}

Real
MoskitoLatHeatIterationXiong::computeDerivative(const Real trail_value, const Real scalar)
{
  Real Uto_plus_tol, Uto_minus_tol, tol_Uto;
  // Derivation concept according to MoskitoEOS2P userobject
  tol_Uto = scalar * _tol;
  Uto_plus_tol = computeResidual(trail_value, scalar + tol_Uto);
  Uto_minus_tol = computeResidual(trail_value, scalar - tol_Uto);

  return (Uto_plus_tol - Uto_minus_tol) / 2.0 / tol_Uto;
}

void
MoskitoLatHeatIterationXiong::computeQpProperties()
{
  if (_rci > _rcem)
    mooseError("Radius casing higher then Radius cement.\n");

  if (_rcem != 0.0 && _rci == 0.0)
    mooseError("Cannot assign Radius cement without Radius casing.\n");

  _RadTubout[_qp] = _rto;
  // Conversion Si to American units
  _TRankine[_qp] = _T[_qp] * gradC_to_gradR;
  _rwb = _diameter[_qp] * m_to_ft / 2.0;

  if (_rwb < _rti)
    mooseError("Wellbore radius smaller then tubing radius.\n");

// Check annulus radii _rao = annulus outer, _rai = annulus inner
if (_rins != 0.0)
  _rai = _rins;
else
  _rai = _rto;

if (_rci != 0.0)
  _rao = _rci;
else if (_rcem != 0.0)
  _rao = _rcem;
else
  _rao = _rwb;

  Real trial = 0.0;
  returnNewtonSolve(_Uto[_qp], _Uto[_qp], _console);
  _Twb[_qp] = TemperatureWBinterface(_Uto[_qp],  _TRankine[_qp]);
}

Real
MoskitoLatHeatIterationXiong::initialGuess(const Real trial_value)
{
  Real ini_guess;
  if (trial_value == 0.0)
    ini_guess = 4.05;
  else
    ini_guess = trial_value;
  return ini_guess;
}
