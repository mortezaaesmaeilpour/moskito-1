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
    type = MoskitoFluidWell1P
    pressure = p
    temperature = 0
    flowrate = q
    well_direction = x
    eos_uo = eos
    viscosity_uo = viscosity
    well_diameter = 0.2
    roughness_type = smooth
    output_properties = 'well_direction_vector density well_velocity well_reynolds_no well_moody_friction'
    gravity = '10 0 0'
    block = 0
  [../]
  [./area1]
    type = MoskitoFluidWell1P
    pressure = p
    temperature = 0
    flowrate = q
    well_direction = x
    eos_uo = eos
    viscosity_uo = viscosity
    well_diameter = 0.25
    roughness_type = smooth
    output_properties = 'well_direction_vector density well_velocity well_reynolds_no well_moody_friction'
    gravity = '10 0 0'
    block = 1
  [../]
[]

[BCs]
  [./pbc]
    type = DirichletBC
    variable = p
    boundary = left
    value = 0
  [../]
  [./qbc]
    type = DirichletBC
    variable = q
    boundary = right
    value = 0
  [../]
[]

[Variables]
  [./p]
  [../]
  [./q]
    scaling = 1e-4
    initial_condition = 0
  [../]
[]

[Kernels]
  [./pkernel]
    type = MoskitoMass1P
    variable = p
    flowrate = q
  [../]
  [./qkernel]
    type = MoskitoMomentum1P
    variable = q
    pressure = p
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
  type = Steady
  l_tol = 1e-10
  l_max_its = 50
  nl_rel_tol = 1e-9
  nl_abs_tol = 1e-10
  nl_max_its = 50
  solve_type = NEWTON
[]

[Postprocessors]
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
[]
