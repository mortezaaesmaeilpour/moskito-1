# "Applied multiphase flow in pipes and flow assurance oil and gas production"
# Al-Safran, E., Brill, J. P., 2017
# Example 2.1: Determine the pipeline inlet pressure?

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
    type = MoskitoPureWater2P
  [../]
  [./viscosity_gas]
    type = MoskitoViscosityConst
    viscosity = 0.00001
  [../]
  [./viscosity_liquid]
    type = MoskitoViscosityConst
    viscosity = 0.001
  [../]
  [./viscosity]
    type = MoskitoViscosity2P
     ve_uo_gas = viscosity_gas
     ve_uo_liquid = viscosity_liquid
  [../]
  [./df]
    type = MoskitoDFHK
  [../]
[]

[Materials]
  [./area0]
    type = MoskitoFluidWell2P
    pressure = p
    enthalpy = h
    flowrate = q
    well_direction = x
    well_diameter = 0.152
    eos_uo = eos
    viscosity_uo = viscosity
    drift_flux_uo = df
    roughness_type = smooth
    gravity = '0 -9.8 0'
    output_properties = 'temperature density well_velocity flow_type_c0 drift_velocity '
    outputs = exodus
  [../]
[]

[BCs]
  [./pbc]
    type = DirichletBC
    variable = p
    boundary = right
    value = '2.07e6'
  [../]
  [./q_right]
      type = MoskitoMassFlowRateCoupled
      boundary = right
      enthalpy = h
      variable = q
      pressure = p
      eos_uo = eos
      mass_flowrate = -46.296296
  [../]
[]

[Variables]
  [./h]
    initial_condition = 83000
  [../]
  [./p]
    initial_condition = 2.07e6
  [../]
  [./q]
    initial_condition = -0.046296296
  [../]
[]

[Kernels]
  [./hkernel]
    type = NullKernel
    variable = h
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
[]
