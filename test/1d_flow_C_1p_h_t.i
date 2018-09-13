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
    flow_rate = q
    well_direction = x
    eos_UO = eos
    viscosity_UO = viscosity
    well_diameter = 0.2
    roughness_type = smooth
    output_properties = 'temperature density well_velocity well_reynolds_no well_moody_friction'
    block = 0
  [../]
  [./area1]
    type = MoskitoFluidWell1P
    pressure = p
    enthalpy = h
    flow_rate = q
    well_direction = x
    eos_UO = eos
    viscosity_UO = viscosity
    well_diameter = 0.35
    roughness_type = smooth
    output_properties = 'temperature density well_velocity well_reynolds_no well_moody_friction'
    block = 1
  [../]
[]

[BCs]
  [./pbc]
    type = DirichletBC
    variable = p
    boundary = left
    value = 2e5
  [../]
  [./qbc]
    type = DirichletBC
    variable = q
    boundary = left
    value = 0.1
  [../]
  [./hbc]
    type = DirichletBC
    variable = h
    boundary = left
    value = 1e5
  [../]
[]

[Variables]
  [./h]
    scaling = 1e-6
    initial_condition = 1e5
  [../]
  [./p]
    # initial_condition = 2e5
  [../]
  [./q]
    scaling = 1e-4
  [../]
[]

[Kernels]
  [./hkernel]
    type = MoskitoEnergy1P
    variable = h
    pressure = p
    flow_rate = q
  [../]
  [./pkernel]
    type = MoskitoMass1P
    variable = p
    flow_rate = q
    enthalpy = h
  [../]
  [./ptkernel]
    type = MoskitoTimeMass1P
    variable = p
    enthalpy = h
  [../]
  [./qkernel]
    type = MoskitoMomentum1P
    variable = q
    pressure = p
    enthalpy = h
  [../]
  [./qtkernel]
    type = MoskitoTimeMomentum1P
    variable = q
    pressure = p
    enthalpy = h
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
  end_time = 20
  num_steps = 10
  l_tol = 1e-10
  l_max_its = 50
  nl_rel_tol = 1e-9
  nl_abs_tol = 1e-10
  nl_max_its = 50
  solve_type = NEWTON
[]

[Postprocessors]
  [./h]
    type = VariableResidual
    variable = h
  [../]
  [./p]
    type = VariableResidual
    variable = p
  [../]
  [./q]
    type = VariableResidual
    variable = q
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
