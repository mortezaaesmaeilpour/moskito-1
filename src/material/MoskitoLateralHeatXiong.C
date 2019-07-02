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

#include "MoskitoLateralHeatXiong.h"

registerMooseObject("MoskitoApp", MoskitoLateralHeatXiong);

template <>
InputParameters
validParams<MoskitoLateralHeatXiong>()
{
    InputParameters params = validParams<Material>();
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
    params.addRequiredParam<Real>("conductivity_earth",
          "Thermal conductivity of the formation (W/(m*K))");
    params.addRequiredParam<Real>("Surface_temperature",
          "Surface temperature (°C); compare with BC");
    params.addRequiredParam<Real>("thermal_diffusivity_earth",
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
          "Raithby & Holland 1975, Churchill 1983");
    MooseEnum time_type
          ("user_time simulation_time");
    params.addParam<MooseEnum>("time_model", time_type,
          "Define time for steady simulations"
          "user_time, simulation_time");
    params.addParam<Real>("user_defined_time", 86400,
          "Time defined by the user for steady state simulation, Default = 1 day");
    return params;
}

MoskitoLateralHeatXiong::MoskitoLateralHeatXiong(const InputParameters & parameters)
  : Material(parameters),
    _T(getMaterialProperty<Real>("temperature")),
    _rti(getParam<Real>("radius_tubbing_inner")),
    _RadTubout(declareProperty<Real>("radius_tubbing_outer")),
    _rto(getParam<Real>("radius_tubbing_outer")),
    _rins(getParam<Real>("radius_insulation")),
    _rci(getParam<Real>("radius_casing_inner")),
    _rcem(getParam<Real>("radius_cement")),
    _eao(getParam<Real>("emissivity_annulus_outer")),
    _eai(getParam<Real>("emissivity_annulus_inner")),
    _rhoA(getParam<Real>("density_annulus")),
    _nuA(getParam<Real>("dyn_viscosity_annulus")),
    _lambdaA(getParam<Real>("conductivity_annulus")),
    _cpA(getParam<Real>("capacity_annulus")),
    _betaA(getParam<Real>("thermal_expansion_annulus")),
    _alphaA(getParam<Real>("thermal_diffusivity_annulus")),
    _lambdaE(getParam<Real>("conductivity_earth")),
    _Twb(declareProperty<Real>("temperature_formation_well_boundary")),
    _alphaE(getParam<Real>("thermal_diffusivity_earth")),
    _Uto(declareProperty<Real>("thermal_resistivity_well")),
    _otU(declareProperty<Real>("lamreht_ytivitsiser_llew")),
    Tsurf(getParam<Real>("Surface_temperature")),
    _diameter(getMaterialProperty<Real>("well_diameter")),
    _lambdaCem(getParam<Real>("conductivity_cement")),
    _lambdaCas(getParam<Real>("conductivity_casing")),
    _lambdaIns(getParam<Real>("conductivity_insulation")),
    _lambdaTub(getParam<Real>("conductivity_tubing")),
    _gravity(getMaterialProperty<RealVectorValue>("gravity")),
    gradT(getFunction("geothermal_gradient")),
    _hc(getParam<MooseEnum>("hc_calucation_model")),
    _td(getParam<MooseEnum>("time_model")),
    _ut(getParam<Real>("user_defined_time"))
{
  Tsurf *= gradC_to_gradR;
  Tsurf += Rankine_absol;
  Boltz = Boltz * Watt_to_Btu_per_h / (m2_to_ft2 * gradC_to_gradR * gradC_to_gradR * gradC_to_gradR * gradC_to_gradR);
  _rti *= m_to_ft;
  _rto *= m_to_ft;
  _rins *= m_to_ft;
  _rci *= m_to_ft;
  _rcem *= m_to_ft;
  _lambdaE *= Watt_to_Btu_per_h / (m_to_ft * gradC_to_gradR);
  _alphaE *= m2_to_ft2;
  _rhoA *= kg_to_lbm / m3_to_ft3;
  _nuA *=  kg_to_lbm / m_to_ft * s_to_h;
  _lambdaA *=  Watt_to_Btu_per_h / (m_to_ft * gradC_to_gradR);
  _cpA *=  J_to_Btu / (kg_to_lbm * gradC_to_gradR);
  _betaA /= gradC_to_gradR;
  _alphaA *= m2_to_ft2 * s_to_h;
  _lambdaCem *=  Watt_to_Btu_per_h / (m_to_ft * gradC_to_gradR);
  _lambdaCas *=  Watt_to_Btu_per_h / (m_to_ft * gradC_to_gradR);
  _lambdaIns *=  Watt_to_Btu_per_h / (m_to_ft * gradC_to_gradR);
  _lambdaTub *=  Watt_to_Btu_per_h / (m_to_ft * gradC_to_gradR);
}



// Compute dimensionless time (ft)
Real
MoskitoLateralHeatXiong::Cal_ft(Real _alphaE, Real _rti, Real _lambdaE, Real _Uto, Real _rto)
{
  Real Timing;
  if (_td == user_time)
    Timing = _ut;
  else
    Timing = _t;

    // std::cout<<"_td = "<<_td<<std::endl;
    // std::cout<<"_ut = "<<_ut<<std::endl;
    // std::cout<<"_t = "<<_t<<std::endl;
    // std::cout<<"Timing = "<<Timing<<std::endl;

  if (Timing < 1209600)
  {
  Real column;
  if(_Uto == 0.0) // Case can only occur in the initial step
    column = _rto / _lambdaE;
  else
    column = _rto * _Uto / _lambdaE;

  Real row = _alphaE * Timing / (_rti * _rti);
  // std::cout<<"column = "<<column<<std::endl;
  // std::cout<<"row = "<<row<<std::endl;
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
    ft = std::log(2 * std::pow(_alphaE * Timing,0.5) / _rwb) - 0.29;

  return ft;
}


Real
MoskitoLateralHeatXiong::Cal_Te(Real Tsurf)
{
  return Tsurf +  gradT.value(_t, _q_point[_qp]) * gradC_to_gradR;
}

void
MoskitoLateralHeatXiong::computeQpProperties()
{
  if (_rci > _rcem)
    mooseError("Radius casing higher then Radius cement.\n");

  if (_rcem != 0.0 && _rci == 0.0)
    mooseError("Cannot assign Radius cement without Radius casing.\n");

  _RadTubout[_qp] = _rto;
  // Conversion Si to American units
  Real _TRankine;
  _TRankine = _T[_qp] * gradC_to_gradR;
  _rwb = _diameter[_qp] * m_to_ft / 2.0;

  if (_rwb < _rti)
    mooseError("Wellbore radius smaller then tubing radius.\n");

// Check annulus radii rao = annulus outer, rai = annulus inner
Real rai, rao;
if (_rins != 0.0)
  rai = _rins;
else
  rai = _rto;

if (_rci != 0.0)
  rao = _rci;
else if (_rcem != 0.0)
  rao = _rcem;
else
  rao = _rwb;

Real ft =  Cal_ft(_alphaE, _rti, _lambdaE, _Uto[_qp], _rto);
Real Test = Cal_Te(Tsurf);
Real Test1 = Tsurf;

// Calculate temperature at cement/formation boundary
_Twb[_qp] = _rto * _Uto[_qp] * ft * _TRankine;
_Twb[_qp] += _lambdaE * Cal_Te(Tsurf);
_Twb[_qp] /= _rto * _Uto[_qp] * ft + _lambdaE;

Real HansPeter = gradT.value(_t, _q_point[_qp]);
// std::cout<<"z coord ="<<_q_point[_qp] <<std::endl;
// std::cout<<"_rto "<<std::setprecision(10)<<_rto<<std::endl;
// std::cout<<"_TRankine "<<std::setprecision(10)<<_TRankine<<std::endl;
// std::cout<<"ft = "<<std::setprecision(10)<<ft<<std::endl;
// std::cout<<"_Uto = "<<std::setprecision(10)<<_Uto[_qp]<<std::endl;
// std::cout<<"_lambdaE = "<<std::setprecision(10)<<_lambdaE<<std::endl;
// std::cout<<" Test = "<<std::setprecision(10)<<Test<<std::endl;
// std::cout<<" Test1 = "<<std::setprecision(10)<<Test1<<std::endl;
// std::cout<<"GRadT "<<std::setprecision(10)<<HansPeter<<std::endl;



// Calculate casing internal temperature
Real Tci = 0.0;
if (_rcem != 0.0)
{
  Tci += std::log(_rwb / _rcem) / _lambdaCem;
  Tci += std::log(_rcem / _rci) / _lambdaCas;
  Tci *= _rto * _Uto[_qp] * (_TRankine - _Twb[_qp]);
}
Tci += _Twb[_qp];

// Calculate tubing external temperature
Real Tto;
Tto = std::log(_rto / _rti) / _lambdaTub;
Tto *= - _rto * _Uto[_qp] * (_TRankine - _Twb[_qp]);
Tto += _TRankine;

// Calculation of Annuls effect - hc
// Calculate Prandtl number
Real Pr;
Pr = _cpA * _nuA / _lambdaA;

// Calculate Grashof number
Real Gr;
Gr = (rao - rai) * (rao - rai) * (rao - rai);
Gr *= _rhoA * _rhoA * _betaA;
Gr *= Tto - Tci;
Gr *= _gravity[_qp].norm() * m_to_ft * s_to_h * s_to_h;
Gr /= _nuA * _nuA;

// Calculate Rayleigh _p_var_number
Real Lc;
Lc = 2 * std::pow(std::log(rao / rai),4.0 / 3.0);
Lc /= std::pow(std::pow(rai,-3.0 / 5.0) + std::pow(rao,-3.0 / 5.0), 5.0 / 3.0);

Real Ray;
Ray = _gravity[_qp].norm() * m_to_ft * s_to_h * s_to_h;
Ray *= _betaA * (Tto - Tci) * Lc * Lc * Lc;
Ray /= _alphaA *  _nuA / _rhoA;

Real khc;
Real hc;
Real fPR;
Real Nu;
switch(_hc)
{
case HC_case::Dropkin_Sommerscales:
// Calculate equivalent thermal conductivity of the annulus fluid
khc = 0.049 * std::pow(Gr * Pr, 0.333);
khc *= std::pow(Pr,0.074) * _lambdaA;
// Calculate convective heat transfer coefficient
hc = khc;
hc /= rai * std::log(rao / rai);
break;

case HC_case::Raithby_Hollands:
// Calculate equivalent thermal conductivity of the annulus fluid
//  I am quite unclear about 0.384. In Xiong it is 0.386 but checking the original document it should be o.48 * 0.8
khc = 0.384 * std::pow(Pr / (0.861 + Pr),0.25);
khc *= std::pow(Ray,0.25) * _lambdaA;
// Calculate convective heat transfer coefficient
hc = khc;
hc /= rai * std::log(rao / rai);
break;


// Document not available
case HC_case::Churchill:
fPR = std::pow(1 + std::pow(0.5 / Pr,9.0 / 16.0), - 16.0 / 9.0);
// Calculate Nusselt number
Nu = 0.364 * std::pow(Ray * fPR, 0.25);
Nu *= std::pow(rao / rai,0.5);
// Calculate convective heat transfer coefficient
hc = Nu * _lambdaA / (2 * rao);
break;
}

// Calculation of Annuls effect - hr
//Calculate auxilary variable
Real OverFtci = 0.0;
OverFtci += 1.0 / _eao - 1.0;
OverFtci *= rai / rao;
OverFtci += 1.0 / _eai;

// Calculate radial heat transfer coefficient
Real hr;
hr = Boltz / OverFtci;
hr *= Tto + Tci;
hr *= Tto * Tto + Tci * Tci;

// Calculate 1 / Uto
_otU[_qp] = std::log(_rto / _rti) / _lambdaTub;
_otU[_qp] += 1.0 / (rai * (hc + hr));
if (_rins != 0.0)
  _otU[_qp] += std::log(_rins / _rto) / _lambdaIns;
if (_rci != 0.0)
  _otU[_qp] += std::log(_rcem / _rci) / _lambdaCas;
if (_rcem != 0.0)
  _otU[_qp] += std::log(_rwb / _rcem) / _lambdaCem;
_otU[_qp] *= _rto;

// Calculate Uto
_Uto[_qp] = 1.0 / _otU[_qp];


// std::cout<<"_TWB = "<<std::setprecision(10)<<_Twb[_qp]<<std::endl;
// std::cout<<"_Uto = "<<std::setprecision(10)<<_Uto[_qp]<<std::endl;
// Real Test = Tto - Tci;
// std::cout<<"Test = "<<std::setprecision(10)<<Test<<std::endl;

// std::cout<<"rao = "<<rao<<std::endl;
// std::cout<<"rai = "<<rai<<std::endl;
// std::cout<<"_rto = "<<_rto<<std::endl;
// std::cout<<"_rti = "<<_rti<<std::endl;
// std::cout<<"_rins = "<<_rins<<std::endl;
// std::cout<<"_rci = "<<_rci<<std::endl;
// std::cout<<"_rcem = "<<_rcem<<std::endl;
// std::cout<<"_rwb = "<<_rwb<<std::endl;
// std::cout<<"_T = "<<_T[_qp]<<std::endl;
// std::cout<<"Tsurf = "<<Tsurf<<std::endl;

// std::cout<<"_lambdaE = "<<_lambdaE<<std::endl;
// std::cout<<"_alphaE = "<<_alphaE<<std::endl;
// std::cout<<"_rhoA = "<<_rhoA<<std::endl;
// std::cout<<"_alphaA = "<<_alphaA<<std::endl;
// std::cout<<"_betaA = "<<_betaA<<std::endl;
// std::cout<<"_nuA = "<<_nuA<<std::endl;
// std::cout<<"_cpA = "<<_cpA<<std::endl;
// std::cout<<"_lambdaA = "<<_lambdaA<<std::endl;
// //
// std::cout<<"_lambdaTub = "<<_lambdaTub<<std::endl;
// std::cout<<"_lambdaCas = "<<_lambdaCas<<std::endl;
// std::cout<<"_lambdaCem = "<<_lambdaCem<<std::endl;
// std::cout<<"_lambdaIns = "<<_lambdaIns<<std::endl;
// std::cout<<"_gradT = "<<Cal_Te(Tsurf)<<std::endl;
//
// std::cout<<"_TRankine = "<<_TRankine<<std::endl;
// std::cout<<"_Twb = "<<std::setprecision(10)<<_Twb[_qp]<<std::endl;
// std::cout<<"Tci = "<<std::setprecision(10)<<Tci<<std::endl;
// std::cout<<"Tto = "<<std::setprecision(10)<<Tto<<std::endl;
// std::cout<<"_nuA = "<<std::setprecision(10)<<_nuA<<std::endl;
// std::cout<<"_rhoA = "<<std::setprecision(10)<<_rhoA<<std::endl;
// std::cout<<"_betaA = "<<std::setprecision(10)<<_betaA<<std::endl;
// std::cout<<"Pr = "<<std::setprecision(10)<<Pr<<std::endl;
// std::cout<<"_gravity[_qp].norm() = "<<std::setprecision(10)<<_gravity[_qp].norm() * m_to_ft * s_to_h * s_to_h<<std::endl;
// std::cout<<"Gr = "<<std::setprecision(10)<<Gr<<std::endl;
// std::cout<<"Lc = "<<std::setprecision(10)<<Lc<<std::endl;
// std::cout<<"Ray = "<<std::setprecision(10)<<Ray<<std::endl;
// std::cout<<"khc = "<<std::setprecision(10)<<khc<<std::endl;
// std::cout<<"hc = "<<std::setprecision(10)<<hc<<std::endl;
// std::cout<<"OverFtci = "<<std::setprecision(10)<<OverFtci<<std::endl;
// std::cout<<"hr = "<<std::setprecision(10)<<hr<<std::endl;
// std::cout<<"fPR = "<<std::setprecision(10)<<fPR<<std::endl;
// std::cout<<"Nu = "<<std::setprecision(10)<<Nu<<std::endl;
// std::cout<<"z coord ="<<_q_point[_qp] <<std::endl;
// std::cout<<"Boltz ="<<Boltz <<std::endl;

// std::cout<<"OverFtci = "<<std::setprecision(10)<<OverFtci<<std::endl;
// std::cout<<"_otU = "<<std::setprecision(10)<<_otU[_qp]<<std::endl;
// std::cout<<"_Uto = "<<std::setprecision(10)<<_Uto[_qp]<<std::endl;


}

// Real angle = std::fabs(_gravity * _well_dir / _gravity.norm());
