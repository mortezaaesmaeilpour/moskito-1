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
    params.addRequiredParam<Real>("radius_wellbore",
          "Radius of the complete wellbore (m)");
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
          "csv file containing the Temperatures at depth z vs T");
    MooseEnum hc_model
          ("Dropkin_Sommerscales Raithby_Hollands Churchill");
    params.addParam<MooseEnum>("hc_calucation_model", hc_model,
          "Type of calculation of hc; Methodology: Dropkin & Sommerscales 1965,"
          "Raithby & Holland 1975, Churchill 1983; While DS and RH are similar in the results, Churchill approach differs strongly");
    params.addParam<Real>("user_defined_time",-1.0,
          "Time defined by the user for steady state simulation, Default = 1 day");
    params.addParam<Real>("derivative_tolerance", 0.001,
          "Tolerance to calculate derivatives based on numerical differentiation");
    params.addParam<RealVectorValue>("Ind_grav", RealVectorValue(9.8,0.0,0.0),
          "Independent gravity to calculate Grashof or Rayleigh number in case overall gravity is set to zero");
    MooseEnum dt_model
          ("Ramey_1962 Kutasov_2003_eq15 Kutasov_2003_eq16 Kutasov_2003_eq17 Kutasov_2003_eq18 Kutasov_1987_eq19 Kutun_2015_eq20");
    params.addParam<MooseEnum>("DimTime_calculation_model", dt_model,
                "Define calculation type of dimensionless time: "
                "Ramey1962, Kutun2015");
    return params;
}

MoskitoLatHeatIterationXiong::MoskitoLatHeatIterationXiong(const InputParameters & parameters)
  : Material(parameters),
    NewtonIteration(parameters),
    _T(getMaterialProperty<Real>("temperature")),
    _rwb(getParam<Real>("radius_wellbore")),
    _RadTubout(declareProperty<Real>("radius_tubbing_outer")),
    _TRock(declareProperty<Real>("formation_temperature")),
    _rto(getParam<Real>("radius_tubbing_outer")),
    _rins(getParam<Real>("radius_insulation")),
    _rci(getParam<Real>("radius_casing_inner")),
    _rcem(getParam<Real>("radius_cement")),
    _Annulus_eao(getParam<Real>("emissivity_annulus_outer")),
    _Annulus_eai(getParam<Real>("emissivity_annulus_inner")),
    _Annulus_rho(getParam<Real>("density_annulus")),
    _Annulus_nu(getParam<Real>("dyn_viscosity_annulus")),
    _Annulus_lambda(getParam<Real>("conductivity_annulus")),
    _Annulus_cp(getParam<Real>("capacity_annulus")),
    _Annulus_beta(getParam<Real>("thermal_expansion_annulus")),
    _Annulus_alpha(getParam<Real>("thermal_diffusivity_annulus")),
    _Rock_lambda(getParam<Real>("conductivity_rock")),
    _Rock_alpha(getParam<Real>("thermal_diffusivity_rock")),
    _Uto(declareProperty<Real>("thermal_resistivity_well")),
    _Twb(declareProperty<Real>("temperature_well_formation_interface")),
    _diameter(getMaterialProperty<Real>("well_diameter")),
    _lambdaCem(getParam<Real>("conductivity_cement")),
    _lambdaCas(getParam<Real>("conductivity_casing")),
    _lambdaIns(getParam<Real>("conductivity_insulation")),
    _lambdaTub(getParam<Real>("conductivity_tubing")),
    _gravity(getMaterialProperty<RealVectorValue>("gravity")),
    _gradT(getFunction("geothermal_gradient")),
    _hc(getParam<MooseEnum>("hc_calucation_model")),
    _ut(getParam<Real>("user_defined_time")),
    _Dim_time(getParam<MooseEnum>("DimTime_calculation_model")),
    _tol(getParam<Real>("derivative_tolerance")),
    _independ_gravity(getParam<RealVectorValue>("Ind_grav")),
    _well_dir(getMaterialProperty<RealVectorValue>("well_direction_vector"))
{
}

// Compute dimensionless time (ft)
Real
MoskitoLatHeatIterationXiong::transienttimefunction(Real Uto)
{
  Real Timing = (_ut == -1.0) ? _t : _ut;

  switch(_Dim_time)
  {
  case Dim_time_case::Ramey_1962:

  if (Timing < 604800)
  {
    Real column = (Uto == 0.0) ? _rto / _Rock_lambda : _rto * Uto / _Rock_lambda;

    Real row = _Rock_alpha * Timing / (_rwb * _rwb);

    int col_number = 0, row_number = 0;

    for(int i = 0; i <= 12; ++i)
    {
      if (column > column_vec[13])
      {
        col_number = 13;
      break;
      }
      else if (column > column_vec[i])
        col_number = i;
      else
      break;
    }

    for(int j = 0; j <= 8; ++j)
    {
      if (row > row_vec[9])
      {
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
    else
    {
      Hor1 = (std::log(column) -  std::log(column_vec[col_number]));
      Hor1 /= (std::log(column_vec[col_number+1]) - std::log(column_vec[col_number]));
      Hor1 *= (Ramey[row_number][col_number+1] - Ramey[row_number][col_number]);
      Hor1 += Ramey[row_number][col_number];
    }

    // Horizontal interpolation 2
    Real Hor2;
    if (col_number == 13)
      Hor2 = Ramey[row_number+1][col_number];
    else
    {
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
    else
    {
      ft = (std::log(row) -  std::log(row_vec[row_number]));
      ft /= (std::log(row_vec[row_number+1]) - std::log(row_vec[row_number]));
      ft *= Hor2 - Hor1;
      ft += Hor1;
    }
  }

  else
    ft = std::log(2.0 * std::pow(_Rock_alpha * Timing,0.5) / _rwb) - 0.29;
  break;

  case Dim_time_case::Kutasov_2003_eq15:
    ft = std::log(1.0 + 1.4986 * std::pow(_Rock_alpha * Timing / (_rwb *_rwb),0.5));
  break;

  case Dim_time_case::Kutasov_2003_eq16:
    ft = std::log(1.0 + (1.781 - 1.0 / (2.701 + std::pow(_Rock_alpha * Timing / (_rwb *_rwb),0.5)) * std::pow(_Rock_alpha * Timing / (_rwb *_rwb),0.5)));
  break;

  case Dim_time_case::Kutasov_2003_eq17:
    ft = 2.0 * std::pow(_Rock_alpha * Timing / (_rwb *_rwb) / PI,0.5);
  break;

  case Dim_time_case::Kutasov_2003_eq18:
    ft = std::log(1.0 + 1.128 * std::pow(_Rock_alpha * Timing / (_rwb *_rwb),0.5));
  break;

  case Dim_time_case::Kutasov_1987_eq19:
    ft = std::log(1.0 + (1.571 - 1.0 / (4.959 + std::pow(_Rock_alpha * Timing / (_rwb *_rwb),0.5)) * std::pow(_Rock_alpha * Timing / (_rwb *_rwb),0.5)));
  break;

  case Dim_time_case::Kutun_2015_eq20:
    ft = std::log(1.0 + 1.7 * std::pow(_Rock_alpha * Timing / (_rwb *_rwb),0.5));
  break;

  }

  return ft;
}

Real
MoskitoLatHeatIterationXiong::TemperatureFormation()
{
  Real Trock = _gradT.value(_t, _q_point[_qp]);
  _TRock[_qp] = Trock;

  return Trock;
}

Real
MoskitoLatHeatIterationXiong::TemperatureWBinterface(Real Uto, Real Temp)
{
  Real Twb = 0.0;
  Twb += _rto * Uto * transienttimefunction(Uto) * Temp;
  Twb += _Rock_lambda * TemperatureFormation();
  Twb /= _rto * Uto * transienttimefunction(Uto) + _Rock_lambda;

  return Twb;
}

Real
MoskitoLatHeatIterationXiong::TemperatureCasingAnnulusInterface(Real Uto, Real Temp)
{
  Real Tci = 0.0;
  if (_rcem != 0.0)
  {
    Tci += std::log(_rwb / _rcem) / _lambdaCem;
    Tci += std::log(_rcem / _rci) / _lambdaCas;
    Tci *= _rto * Uto * (Temp - TemperatureWBinterface(Uto, Temp));
  }
  Tci += TemperatureWBinterface(Uto, Temp);

  return Tci;
}

Real
MoskitoLatHeatIterationXiong::TemperatureTubingOuter(Real Uto, Real Temp)
{
  Real Tto;
  Tto = std::log(_rto / _rti) / _lambdaTub;
  Tto *= - _rto * Uto * (Temp - TemperatureWBinterface(Uto, Temp));
  Tto += Temp;

  return Tto;
}

Real
MoskitoLatHeatIterationXiong::RadialHeatTransferCoefficient(Real Uto, Real Temp)
{
  Real OverFtci = 0.0;
  OverFtci += 1.0 / _Annulus_eao - 1.0;
  OverFtci *= _rto / _rci;
  OverFtci += 1.0 / _Annulus_eai;

  Real hr;
  hr = Boltz / OverFtci;
  hr *= TemperatureTubingOuter(Uto, Temp) + TemperatureCasingAnnulusInterface(Uto, Temp);
  hr *= TemperatureTubingOuter(Uto, Temp) * TemperatureTubingOuter(Uto, Temp) + TemperatureCasingAnnulusInterface(Uto, Temp) * TemperatureCasingAnnulusInterface(Uto, Temp);

  return hr;
}

Real
MoskitoLatHeatIterationXiong::ConvectiveHeatTransferCoefficient(Real Uto, Real Temp, Real grav)
{
  Real khc, hc;
  // Calculate Prandtl number
  Real Pr;
  Pr = _Annulus_cp * _Annulus_nu / _Annulus_lambda;

  switch(_hc)
  {
  case HC_case::Dropkin_Sommerscales:

  khc = 0.049 * std::pow(Grashof(Uto, Temp, grav) * Pr, 0.333);
  khc *= std::pow(Pr,0.074) * _Annulus_lambda;

  hc = khc / (_rai * std::log(_rao / _rai));
  break;

  case HC_case::Raithby_Hollands:
  //  I am quite unclear about 0.384. In Xiong it is 0.386 but checking the original document it should be 0.48 * 0.8
  khc = 0.384 * std::pow(Pr / (0.861 + Pr),0.25);
  khc *= std::pow(Rayleigh(grav, Uto, Temp),0.25) * _Annulus_lambda;

  hc = khc / (_rai * std::log(_rao / _rai));
  break;

  // Document not available
  case HC_case::Churchill:

  Real fPR, Nu;
  fPR = std::pow(1 + std::pow(0.5 / Pr,9.0 / 16.0), - 16.0 / 9.0);
  // Calculate Nusselt number
  Nu = 0.364 * std::pow(Rayleigh(grav, Uto, Temp) * fPR, 0.25);
  Nu *= std::pow(_rao / _rai,0.5);

  hc = Nu * _Annulus_lambda / (2.0 * _rao);
  break;
  }

  return hc;
}

Real
MoskitoLatHeatIterationXiong::Grashof(Real Uto, Real Temp, Real grav)
{
  Real Gr;

  Gr = (_rao - _rai) * (_rao - _rai) * (_rao - _rai);
  Gr *= _Annulus_rho * _Annulus_rho * _Annulus_beta;
  Gr *= std::abs(TemperatureTubingOuter(Uto, Temp) - TemperatureCasingAnnulusInterface(Uto, Temp));
  Gr *= grav;
  Gr /= _Annulus_nu * _Annulus_nu;

  return Gr;
}

Real
MoskitoLatHeatIterationXiong::Rayleigh(Real grav, Real Uto, Real Temp)
{
  Real Lc, Ray;
  Lc = 2 * std::pow(std::log(_rao / _rai),4.0 / 3.0);
  Lc /= std::pow(std::pow(_rai,-3.0 / 5.0) + std::pow(_rao,-3.0 / 5.0), 5.0 / 3.0);
  Ray = grav * Lc * Lc * Lc;
  Ray *= _Annulus_beta * std::abs(TemperatureTubingOuter(Uto, Temp) - TemperatureCasingAnnulusInterface(Uto, Temp));
  Ray /= _Annulus_alpha *  _Annulus_nu / _Annulus_rho;

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
  hr =  RadialHeatTransferCoefficient(scalar, _T[_qp]);
  Real hc, grav;
  grav = _gravity[_qp] * _well_dir[_qp];

  if (grav == 0.0)
    grav = _independ_gravity * _well_dir[_qp];

  hc = ConvectiveHeatTransferCoefficient(scalar, _T[_qp], grav);

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
  _rti = _diameter[_qp] / 2.0;

  if (_rwb < _rti)
    mooseError("Wellbore radius smaller then tubing radius.\n");

// Check annulus radii _rao = annulus outer, _rai = annulus inner
_rai = (_rins != 0.0) ? _rins : _rto;

if (_rci != 0.0)
  _rao = _rci;
else if (_rcem != 0.0)
  _rao = _rcem;
else
  _rao = _rwb;

  Real trial = 0.0;
  returnNewtonSolve(_Uto[_qp], _Uto[_qp], _console);
  _Twb[_qp] = TemperatureWBinterface(_Uto[_qp],  _T[_qp]);
}

Real
MoskitoLatHeatIterationXiong::initialGuess(const Real trial_value)
{
  Real ini_guess;
  if (trial_value == 0.0)
    ini_guess = 22.99696519;
  else
    ini_guess = trial_value;
  return ini_guess;
}
