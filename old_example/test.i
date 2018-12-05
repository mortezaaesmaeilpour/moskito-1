[Mesh]
  type = GeneratedMesh
  dim = 1
  xmin = 0
  xmax = 1
  nx = 1000
[]

[UserObjects]
  [./eos_g]
    type = MoskitoEOSIdealGas
    molar_mass = 2.16e-2
  [../]
  [./eos_l]
    type = MoskitoEOSIdealFluid
  [../]
  [./eos_2p]
    type = MoskitoEOS2P
    eos_uo_gas = eos_g
    eos_uo_liquid = eos_l
  [../]
  [./ve_g]
    type = MoskitoViscosityConst
    viscosity = 0.00001
  [../]
  [./ve_l]
    type = MoskitoViscosityConst
    viscosity = 0.001
  [../]
  [./ve_2p]
    type = MoskitoViscosity2P
    ve_uo_gas = ve_g
    ve_uo_liquid = ve_l
    mixing_type = Mean_ME12
  [../]
[]

[Functions]
  [./void]
    type = ParsedFunction
    value = 'x'
  [../]
[]

[Materials]
  [./g]
    type = GenericConstantMaterial
    prop_names = 'drift_velocity flow_type_c0'
    prop_values = '0.0 1.0'
  [../]
  [./v]
    type = GenericFunctionMaterial
    prop_names = 'void_fraction'
    prop_values = 'void'
    outputs = exodus
  [../]
  [./area0]
    type = MoskitoFluidWell2P
    pressure = p
    enthalpy = h
    flowrate = q
    well_direction = x
    eos_uo = eos_2p
    viscosity_uo = ve_2p
    well_diameter = 0.152
    roughness_type = smooth
    gravity = '0 -9.8 0'
    outputs = exodus
  [../]
[]

[Variables]
  [./h]
  [../]
  [./p]
  [../]
  [./q]
  [../]
[]

[Kernels]
  [./hkernel]
    type = NullKernel
    variable = h
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
