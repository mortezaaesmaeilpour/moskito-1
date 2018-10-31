[Mesh]
  type = GeneratedMesh
  dim = 1
  xmin = 0
  xmax = 1000
  nx = 1000
[]

[MeshModifiers]
  [./block1]
    type = SubdomainBoundingBox
    block_id = 1
    bottom_left = '-1 -1 -1'
    top_right = '500 1 1'
  [../]
[]

[UserObjects]
  [./eos]
    type = MoskitoEOSIdealFluid
    reference_pressure = 0
    thermal_conductivity = 10000
  [../]
  [./viscosity]
    type = MoskitoViscosityConst
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
    output_properties = 'temperature density well_velocity well_reynolds_no well_moody_friction'
    block = 0
  [../]
  [./area1]
    type = MoskitoFluidWell1P
    pressure = p
    enthalpy = h
    flowrate = q
    well_direction = x
    eos_uo = eos
    viscosity_uo = viscosity
    well_diameter = 0.35
    roughness_type = smooth
    output_properties = 'temperature density well_velocity well_reynolds_no well_moody_friction'
    block = 1
  [../]
[]

[BCs]
  [./hbc]
    type = MoskitoEnthalpyTemperatureDBC
    variable = h
    boundary = left
    temperature = 300
    eos_uo = eos
  [../]
[]

[Variables]
  [./h]
    # scaling = 1e-6
    initial_condition = 83950
  [../]
  [./p]
    initial_condition = 0
  # initial_condition = 2e5
  [../]
  [./q]
    initial_condition = 0
    # scaling = 1e-6
  [../]
[]

[Kernels]
  [./hkernel]
    type = MoskitoEnergy1P
    variable = h
    pressure = p
    flowrate = q
  [../]
  [./htkernel]
    type = MoskitoTimeEnergy1P
    variable = h
    pressure = p
    flowrate = q
  [../]
  [./pkernel]
    type = NullKernel
    variable = p
  [../]
  [./qkernel]
    type = NullKernel
    variable = q
  [../]
[]

[Preconditioning]
  active = 'p3'
  [./p1]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_hypre_type'
    petsc_options_value = 'hypre boomeramg'
  [../]
  [./p2]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type -ksp_gmres_restart'
    petsc_options_value = 'asm lu NONZERO 51'
  [../]
  [./p3]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -ksp_type -sub_pc_type -pc_asm_overlap -sub_pc_factor_shift_type -ksp_gmres_restart'
    petsc_options_value = 'asm gmres lu 2 NONZERO 51'
  [../]
[]

[Executioner]
  type = Transient
  end_time = 10000000
  # num_steps = 50
  l_tol = 1e-8
  l_max_its = 50
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  nl_max_its = 50
  solve_type = NEWTON
  [./TimeStepper]
    type = SolutionTimeAdaptiveDT
    dt = 1000
  [../]
[]

[Outputs]
  # exodus = true
  print_linear_residuals = true
  [./test]
    type = Exodus
    output_material_properties = true
  [../]
  # [./debug]
  # type = VariableResidualNormsDebugOutput
  # output_nonlinear = true
  # [../]
[]
