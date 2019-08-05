
[Mesh]
  type = GeneratedMesh
  dim = 1
  xmin = 100
  xmax = 1260
  nx = 2
[]

# [MeshModifiers]
#   [./rotate]
#     type = Transform
#     transform = ROTATE
#     vector_value = '7 0 0'
#   [../]
# []

[UserObjects]
  [./eos]
    type = MoskitoEOSIdealFluid
    bulk_modulus = 2e+012
  [../]
  [./viscosity]
    type = MoskitoViscosityConst
    viscosity = 0.001
  [../]
[]

[Functions]
  [./grad_func]
    type = ParsedFunction
    value = '0.03 * x'
  [../]
[]

[Materials]
  [./area0]
    type = MoskitoFluidWell1P
    pressure = p
    enthalpy = h
    flowrate = q
    well_direction = x
    eos_uo = eos
    viscosity_uo = viscosity
    well_diameter = 0.2
    roughness_type = smooth
    gravity = '9.8 0 0'
    output_properties = 'temperature density well_velocity specific_heat'
    outputs = exodus
  [../]
  [./LateralHeat]
    type = MoskitoLateralHeatXiong
    internal_solve_output_on = always
    radius_tubbing_inner = 0.051
    radius_tubbing_outer = 0.055
    # radius_casing_inner = 0.108204003
    # radius_cement = 0.121920004
    Surface_temperature = 15
    conductivity_earth = 2.5
    thermal_diffusivity_earth = 1.15e-6
    # conductivity_annulus = 0.6
    # density_annulus = 1000
    # dyn_viscosity_annulus = 0.0000285231
    capacity_annulus = 1025.7659436
    # thermal_expansion_annulus = 0.001755
    conductivity_tubing = 80.42534006
    # conductivity_casing = 80.42534006
    # conductivity_cement = 0.34591544112
    geothermal_gradient = grad_func
    hc_calucation_model = Raithby_Hollands
    emissivity_annulus_outer = 0.9
    emissivity_annulus_inner = 0.9
    user_defined_time = 6666666
    time_model = user_time
    output_properties = 'thermal_resistivity_well hc_calucation_model temperature_formation_well_boundary'
    outputs = exodus
  [../]
[]

[BCs]
  [./pbc]
    type = DirichletBC
    variable = p
    boundary = right
    value = '5.81e6'
  [../]
  [./qbc]
    type = DirichletBC
    variable = q
    boundary = left
    value = 0.01708
  [../]
  [./hbc]
    type = DirichletBC
    variable = h
    boundary = right
    value = 1000000
  [../]
[]

[Variables]
  [./h]
    initial_condition = 1000000
    scaling = 1e-3
  [../]
  [./p]
    initial_condition = 5.81e6
     scaling = 1e5
  [../]
  [./q]
    # scaling = 1e-4
    initial_condition = 0.01708
  [../]
[]

[Kernels]
  [./hkernel]
    type = MoskitoEnergy
    variable = h
    pressure = p
    flowrate = q
  [../]
  [./pkernel]
    type = MoskitoMass
    variable = p
    flowrate = q
    enthalpy = h
  [../]
  [./qkernel]
    type = MoskitoMomentum
    variable = q
    pressure = p
    enthalpy = h
  [../]
  [./Heat]
    type = MoskitoHeat
    variable = h
  [../]
[]

[Preconditioning]
  [./p2]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type -ksp_gmres_restart'
    petsc_options_value = 'asm lu NONZERO 51'
  [../]
[]

[Executioner]
  type = Steady
  l_tol = 1e-10
  l_max_its = 50
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-9
  nl_max_its = 50
  solve_type = NEWTON
[]

[Outputs]
  exodus = true
  [./test]
    type = VariableResidualNormsDebugOutput
  [../]
 []
