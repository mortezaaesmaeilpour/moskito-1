[Mesh]
  type = FileMesh
  file = example_2_2.msh
  uniform_refine = 2
[]

[UserObjects]
  [./eos]
    type = MoskitoEOSIdealFluid
    density0 = 883
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
    well_diameter = 0.1016
    roughness_type = smooth
    manual_friction_factor = 0.02
    output_properties = 'well_direction_vector density well_velocity well_reynolds_no well_moody_friction'
  [../]
[]

[BCs]
  [./pbcl]
    type = DirichletBC
    variable = p
    boundary = left
    value = 0
  [../]
  [./qbc]
    type = DirichletBC
    variable = q
    boundary = left
    value = 0.00223
  [../]
[]

[Variables]
  [./p]
  [../]
  [./q]
    scaling = 1e-5
    initial_condition = 0.00223
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
