[Mesh]
  type = GeneratedMesh
  dim = 1
  xmin = 0
  xmax = 248.4
  nx = 100
[]

[MeshModifiers]
  [./openhole]
    type = SubdomainBoundingBox
    block_id = 1
    bottom_left = '159.6 -1 0'
    top_right =   '250 1 0'
  [../]
[]

[UserObjects]
  [./viscosity_gas]
    type = MoskitoViscosityWaterSmith
  [../]
  [./viscosity_liqid]
    type = MoskitoViscosityWaterVogel
  [../]
  [./viscosity_2p]
    type = MoskitoViscosity2P
    ve_uo_gas = viscosity_gas
    ve_uo_liquid = viscosity_liqid
  [../]
  [./df]
    type = MoskitoDFHK
    surface_tension = 0.0288
  [../]
  [./eos]
    type = MoskitoPureWater2P
    derivative_tolerance = 1e-3
    execute_on = LINEAR
  [../]
  [./init_uo]
    type = SolutionUserObject
    mesh = Prep_Garg_inp.e
    timestep = 2
    system_variables = 'p q'
    execute_on = INITIAL
  [../]
[]

[Functions]
  [./pres]
    type = SolutionFunction
    solution = init_uo
    from_variable = 'p'
  [../]
  [./flow]
    type = SolutionFunction
    solution = init_uo
    from_variable = 'q'
  [../]
[]


[GlobalParams]
  pressure = p
  enthalpy = h
  flowrate = q
  well_direction = x
  eos_uo = eos
  viscosity_uo = viscosity_2p
  drift_flux_uo = df
  roughness_type = smooth
  gravity = '9.8 0 0'
  outputs = exodus
  output_properties = 'gas_velocity liquid_velocity void_fraction mass_fraction flow_pattern current_phase gas_density liquid_density density temperature well_velocity'
[]


[Materials]
  [./area]
    type = MoskitoFluidWell2P
    well_diameter = 0.1
  [../]
[]

[BCs]
  [./pbc]
    type = DirichletBC
    variable = p
    boundary = right
    value = 1.918e6
  [../]
  [./qbc]
    type = DirichletBC
    variable = q
    boundary = left
    value = 0.008
  [../]
  [./hbc]
    type = DirichletBC
    variable = h
    boundary = right
    value = 7e5
  [../]
[]

# [AuxVariables]
#   [./h]
#     initial_condition = 7e5
#   [../]
# []

[Variables]
  [./h]
    [./InitialCondition]
      type = FunctionIC
      variable = h
      function = 7e5
    [../]
    scaling = 1e-4
  [../]
  [./p]
    [./InitialCondition]
      type = FunctionIC
      variable = p
      function = pres
    [../]
  [../]
  [./q]
    [./InitialCondition]
      type = FunctionIC
      variable = q
      function = flow
    [../]
    scaling = 1e-6
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
[]

[Preconditioning]
  active = pn1
  [./p1]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type'
    petsc_options_value = ' bjacobi  ilu          NONZERO                 '
  [../]
  [./pn1]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type -snes_type -snes_linesearch_type'
    petsc_options_value = ' bjacobi  ilu          NONZERO                   newtontr   basic               '
  [../]
  # Newton method (no JFNK and PJFNK)
  [./n_p] # (parallel)
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -snes_type -snes_linesearch_type'
    petsc_options_value = ' lu       mumps                         newtonls   basic               '
  [../]
  [./n_s] # Linear timesteps vanished (serial)
    type = SMP
    full = true
    #petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-ksp_type -pc_type -snes_type -snes_linesearch_type'
    petsc_options_value = ' preonly   lu       newtonls   basic               '
  [../]
  # JFNK and PJFNK
  [./k_p1] # Basic (parrallel)
    type = SMP
    full = true
    #petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type -sub_pc_factor_shift_type -ksp_gmres_restart -snes_type -snes_linesearch_type'
    petsc_options_value = ' gmres     hypre    boomeramg      NONZERO                   51                 newtonls   basic               '
  [../]
  [./k_p2] # Middle (serial)
    type = SMP
    full = true
    #petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -pc_asm_overlap -sub_pc_factor_shift_type -ksp_gmres_restart -snes_type -snes_linesearch_type'
    petsc_options_value = ' gmres     asm      lu           2               NONZERO                   51                 newtonls   basic               '
  [../]
[]

[Executioner]
  type = Steady
  l_max_its = 50
  nl_max_its = 50
  l_tol = 1e-8
  nl_rel_tol = 1e-8
  # nl_abs_tol = 1e-1
  # solve_type = NEWTON
  # automatic_scaling = true
[]

[Outputs]
  exodus = true
  [./test]
    type = VariableResidualNormsDebugOutput
  [../]
[]
