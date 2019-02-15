# "Applied multiphase flow in pipes and flow assurance oil and gas production"
# Al-Safran, E., Brill, J. P., 2017
# Example 2.1: Determine the pipeline inlet pressure?

[Mesh]
  type = GeneratedMesh
  dim = 1
  xmin = 0
  xmax = 90
  nx = 9
[]

[UserObjects]
  [./viscosity_gas]
    type = MoskitoViscosityConst
    viscosity = 0.001
  [../]
  [./viscosity_liqid]
    type = MoskitoViscosityConst
    viscosity = 0.01
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
    eos_uo = eos
  [../]
[]

[Variables]
  [./h]
    initial_condition = 10000
  [../]
  [./p]
    initial_condition = 1
  [../]
  [./q]
    initial_condition = 0
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
