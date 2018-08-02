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
[]

[Materials]
  [./area0]
    type = MoskitoSinglePhaseFluidWell
    temperature = 0
    density = rho
    flow_rate = q
    well_direction = x
    eos_UO = eos
    well_diameter = 0.1016
    roughness_type = smooth
    output_properties = 'well_direction_vector pressure_difference well_velocity well_reynolds_no well_moody_friction'
  [../]
[]

[BCs]
  [./rhobcl]
    type = MoskitoDensityCoupledBC
    variable = rho
    boundary = left
    eos_UO = eos
    temperature = 0
    pressure = '0'
  [../]
  [./qbc]
    type = DirichletBC
    variable = q
    boundary = left
    value = 0.00223
  [../]
[]

[Variables]
  [./rho]
    initial_condition = 883
  [../]
  [./q]
    scaling = 1e-7
    initial_condition = 0.00223
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
