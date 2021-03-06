[Mesh]
  type = GeneratedMesh
  dim = 1
  xmin = 0
  xmax = 3000
  nx = 300
[]

[UserObjects]
  [./viscosity_gas]
    type = MoskitoViscosityConst
    viscosity = 0.0001
  [../]
  [./viscosity_liqid]
    type = MoskitoViscosityConst
    viscosity = 0.001
  [../]
  [./viscosity_2p]
    type = MoskitoViscosity2P
    ve_uo_gas = viscosity_gas
    ve_uo_liquid = viscosity_liqid
  [../]
  [./df]
    type = MoskitoDFHK
  [../]
  [./eos]
    type = MoskitoPureWater2P
  [../]
[]

[Materials]
  [./area1]
    type = MoskitoFluidWell2P
    drift_flux_uo = df
    viscosity_uo = viscosity_2p
    pressure = p
    enthalpy = h
    flowrate = q
    well_direction = x
    well_diameter = 0.152
    roughness_type = smooth
    gravity = '9.8 0 0'
    outputs = exodus
    output_properties = 'gas_density liquid_density mass_fraction density specific_heat temperature'
    eos_uo = eos
  [../]
[]

[Variables]
  [./h]
    # initial_condition = 0.7e6
    [./InitialCondition]
      type = FunctionIC
      function = 600000+1000*x
      variable = h
    [../]
  [../]
  [./p]
    initial_condition = 101325
    # [./InitialCondition]
    #   type = FunctionIC
    #   function = 101325+1000*x
    #   variable = p
    # [../]
  [../]
  [./q]
    initial_condition = 0
  [../]
[]

[BCs]
  [./p]
    type = DirichletBC
    variable = p
    boundary = left
    value = 101325
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
    enthalpy = h
    flowrate = q
  [../]
  [./qkernel]
    type = NullKernel
    variable = q
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
