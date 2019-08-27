# Validation of Heat exchange using the Paper from Satman & Tureyen 2016
[Mesh]
  type = GeneratedMesh
  dim = 1
  xmin = 0
  xmax = 2500
  nx = 400
[]

[UserObjects]
  [./viscosity]
    type = MoskitoViscositySmith
  [../]
  [./eos]
    type = MoskitoEOSIdealFluid
    specific_heat = 3160
  [../]
[]

[Functions]
  [./grad_func]
    type = ParsedFunction
    value = '0.09 * x'
  [../]
[]

[Materials]
  [./area0]
    type = MoskitoFluidWell1P
    pressure = p
    enthalpy = h
    flowrate = q
    well_direction = x
    well_diameter = 0.3
    eos_uo = eos
    viscosity_uo = viscosity
    roughness_type = smooth
    gravity = '0.0 0 0'
    outputs = exodus
    output_properties = 'temperature density well_velocity specific_heat well_reynolds_no well_moody_friction viscosity diameter'
  [../]
  [./Lateral]
    type = MoskitoLatHeatIterationXiongSI
     radius_wellbore = 0.150000000002
     radius_tubbing_outer = 0.150000000001
     conductivity_tubing = 40
     # Rock parameters
     Surface_temperature = 293.15
     thermal_diffusivity_rock = 1.102e-6
     conductivity_rock = 2.92
     # Configuration of material
     geothermal_gradient = grad_func
     hc_calucation_model = Dropkin_Sommerscales
     time_model = user_time
     user_defined_time = 8640
     # user_defined_time = 2592000
     DimTime_calculation_model = Kutun_2015
     # internal_solve_full_iteration_history = true
     output_properties = 'Temperature_Rankine temperature_well_formation_interface heat_loss formation_temperature'
     outputs = exodus
   [../]
[]

[Variables]
  [./h]
    scaling = 1e-6
    initial_condition = 84010
  [../]
  [./p]
    initial_condition = 1.0e5
    scaling = 1e1
  [../]
  [./q]
    initial_condition = 0.020038
    scaling = 1e-2
  [../]
[]

[BCs]
  [./qbc]
    type = DirichletBC
    variable = q
    boundary = right
    value = -0.020038
  [../]
  [./pbc]
    type = DirichletBC
    variable = p
    boundary = left
    value = 1.0e5
  [../]
  [./hbc]
    type = MoskitoTemperatureToEnthalpy1P
    variable = h
    pressure = p
    boundary = right
    temperature = 518.15
    eos_uo = eos
  [../]
[../]

[Kernels]
  [./hkernel]
    type = MoskitoEnergy
    variable = h
    pressure = p
    flowrate = q
  [../]
  [./qkernel1]
    type = MoskitoMomentum
    variable = q
    pressure = p
    enthalpy = h
  [../]
  [./pkernel1]
    type = MoskitoMass
    variable = p
    enthalpy = h
    flowrate = q
  [../]
  [./heat]
    type = MoskitoHeat
    variable = h
  [../]
[]

[Preconditioning]
  active = p3
  [./p1]
    type = SMP
    full = true
    #petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-pc_type -pc_hypre_type -snes_type -snes_linesearch_type -sub_pc_factor_shift_type'
    petsc_options_value = 'hypre boomeramg newtonls basic NONZERO'
  [../]
  [./p2]
    type = SMP
    full = true
    #petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-pc_type -sub_pc_type -snes_type -snes_linesearch_type -sub_pc_factor_shift_type -ksp_gmres_restart'
    petsc_options_value = 'asm lu newtonls basic NONZERO 51'
  [../]
  [./p3]
    type = SMP
    full = true
    #petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-pc_type -ksp_type -sub_pc_type -snes_type -snes_linesearch_type -pc_asm_overlap -sub_pc_factor_shift_type -ksp_gmres_restart'
    petsc_options_value = 'asm gmres lu newtonls basic 2 NONZERO 51'
  [../]
[]

[Executioner]
  type = Steady
  l_max_its = 50
  l_tol = 1e-10
  nl_rel_tol = 1e-8
  nl_max_its = 50
  solve_type = NEWTON
  nl_abs_tol = 1e-7
[]

[Outputs]
  exodus = true
  print_linear_residuals = true
  # [./test]
  #   type = VariableResidualNormsDebugOutput
  # [../]
  console = true
[]
