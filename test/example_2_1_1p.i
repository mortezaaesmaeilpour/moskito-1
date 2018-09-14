[Mesh]
  type = GeneratedMesh
  dim = 1
  xmin = 0
  xmax = 9000
  nx = 900
[]

[MeshModifiers]
  [./rotate]
    type = Transform
    transform = ROTATE
    vector_value = '7 0 0'
  [../]
[]

[UserObjects]
  [./eos]
    type = MoskitoEOSIdealFluid
    bulk_modulus = 2e+012
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
    eos_UO = eos
    viscosity_UO = viscosity
    well_diameter = 0.152
    roughness_type = smooth
    output_properties = 'well_direction_vector density well_velocity well_reynolds_no well_moody_friction'
  [../]
[]

[BCs]
  [./pbc]
    type = DirichletBC
    variable = p
    boundary = right
    value = '2.07e6'
  [../]
  [./qbc]
    type = DirichletBC
    variable = q
    boundary = right
    value = 0.046296296
  [../]
[]

[Variables]
  [./p]
    initial_condition = 2.07e6
  [../]
  [./q]
    scaling = 1e-5
    initial_condition = 0.046296296
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
    gravity = '0 -9.8 0'
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
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-9
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
