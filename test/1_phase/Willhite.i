# Over-all heat transfer coefficients in Steam and hot water injection wells
# Willhite, G. P. 19
# Appendix: Sample Calculation - Results of Uto differs slighly, because of rounding of the author

# Example has a lenght that is not of importance
[Mesh]
  type = GeneratedMesh
  dim = 1
  xmin = 0
  xmax = 1
  nx = 1
[]

[UserObjects]
  [./eos]
    type = MoskitoPureWater2P
  [../]
  [./viscosity_gas]
    type = MoskitoViscosityConst
    viscosity = 0.00001
  [../]
  [./viscosity_liquid]
    type = MoskitoViscosityConst
    viscosity = 0.001
  [../]
  [./viscosity]
    type = MoskitoViscosity2P
     ve_uo_gas = viscosity_gas
     ve_uo_liquid = viscosity_liquid
  [../]
  [./df]
    type = MoskitoDFHK
    surface_tension = 1
  [../]
[]

# No temperature gradient used in the example
[Functions]
  [./grad_func]
    type = ParsedFunction
    value = '0.0 * x'
  [../]
[]

[Materials]
  [./area0]
    type = MoskitoFluidWell2P
    pressure = p
    enthalpy = h
    flowrate = q
    well_direction = x
    well_diameter = 0.0890016
    eos_uo = eos
    viscosity_uo = viscosity
    drift_flux_uo = df
    roughness_type = smooth
    gravity = '9.8 0 0'
    outputs = exodus
  [../]
  [./Lateral]
    type = MoskitoLatHeatIterationXiong
     # Geometry of the well. As the example did not contain any tubing radius, which is required for teh material it was artificially set to a small radius
     radius_tubbing_outer = 0.044500801
     radius_casing_inner = 0.108204
     radius_cement = 0.12192
     radius_wellbore = 0.1524
     conductivity_cement = 0.346146923
     conductivity_tubing = 80.42534006
     conductivity_casing = 80.42534006
     # Rock parameters
     Surface_temperature = 37.7778
     thermal_diffusivity_rock = 0.000000738063
     conductivity_rock = 1.7307346
     # Annulus parameters representing a stagnant gas at 14.7psia @ 296°C
     density_annulus = 0.6215162
     conductivity_annulus = 0.0441337
     dyn_viscosity_annulus = 0.0000285262
     capacity_annulus = 1025.76594
     thermal_expansion_annulus = 1.755e-3
     emissivity_annulus_outer = 0.9
     emissivity_annulus_inner = 0.9
     # Configuration of material
     geothermal_gradient = grad_func
     hc_calucation_model = Dropkin_Sommerscales
     user_defined_time = 1814400
     time_model = user_time
     internal_solve_full_iteration_history = true
   [../]
[]

[Variables]
  [./h]
    initial_condition = 3084870
  [../]
  [./p]
    initial_condition = 1.0e6
  [../]
  [./q]
    initial_condition = -0.01
  [../]
[]

[BCs]
  [./hbc]
    type = DirichletBC
    variable = h
    boundary = left
    value = 3084870
  [../]
[../]

[Kernels]
  [./hkernel]
    type = MoskitoEnergy
    variable = h
    pressure = p
    flowrate = q
  [../]
  [./pkernel1]
    type = NullKernel
    variable = p
    flowrate = q
    enthalpy = h
  [../]
  [./qkernel1]
    type = NullKernel
    variable = q
    enthalpy = h
    pressure = p
  [../]
  [./heat]
    type = MoskitoHeat
    variable = h
  [../]
[]

[Preconditioning]
  active = p1
  [./p1]
    type = SMP
    full = true
    #petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-pc_type -pc_hypre_type -snes_type -snes_linesearch_type -sub_pc_factor_shift_type'
    petsc_options_value = 'hypre boomeramg newtonls basic NONZERO'
  [../]
[]

[Executioner]
  type = Steady
  l_max_its = 2
  nl_rel_tol = 1e-4
  nl_abs_tol = 1e-4
  nl_max_its = 50
  solve_type = NEWTON
[]

[Outputs]
  exodus = true
[]
