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
    bulk_modulus = 2e+09
  [../]
  [./viscosity]
    type = MoskitoViscosityConst
  [../]
[]

[Materials]
  [./area0]
    type = MoskitoSinglePhaseFluidWell
    temperature = 0
    density = rho
    flow_rate = q
    well_direction = x
    eos_UO = eos
    viscosity_UO = viscosity
    well_diameter = 0.2
    roughness_type = smooth
    output_properties = 'well_direction_vector pressure_difference well_velocity well_reynolds_no'
    block = 0
  [../]
  [./area1]
    type = MoskitoSinglePhaseFluidWell
    temperature = 0
    well_diameter = 0.25
    flow_rate = q
    roughness_type = smooth
    output_properties = 'well_direction_vector pressure_difference well_velocity well_reynolds_no'
    well_direction = x
    density = rho
    eos_UO = eos
    viscosity_UO = viscosity
    block = 1
  [../]
[]

[BCs]
  [./rhobc]
    type = DirichletBC
    variable = rho
    boundary = left
    value = 1000
  [../]
  [./qbc]
    type = DirichletBC
    variable = q
    boundary = left
    value = 0
  [../]
[]

[Variables]
  [./rho]
    initial_condition = 1000
  [../]
  [./q]
    scaling = 1e-7
    initial_condition = 0
  [../]
[]

[Kernels]
  [./rhokernel]
    type = MoskitoCMass
    variable = rho
    flow_rate = q
  [../]
  [./qkernel]
    type = MoskitoMomentum
    variable = q
    density = rho
    gravity = '10 0 0'
  [../]
[]

[Preconditioning]
  active = 'p2'
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
  type = Steady
  l_tol = 1e-10
  l_max_its = 50
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  nl_max_its = 50
  solve_type = NEWTON
[]

[Postprocessors]
  [./rho]
    type = VariableResidual
    variable = rho
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
[]
