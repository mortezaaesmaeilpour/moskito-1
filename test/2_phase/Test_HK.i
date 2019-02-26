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
    vector_value = '83 0 0'
  [../]
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
    surface_tension = 0.1
  [../]
  [./eos]
    type = MoskitoPureWater2P
  [../]
[]

[Materials]
  [./area0]
    type = MoskitoFluidWell2P
    pressure = p
    enthalpy = h
    flowrate = q
    well_direction = x
    eos_uo = eos
    viscosity_uo = viscosity_2p
    drift_flux_uo = df
    well_diameter = 0.152
    roughness_type = smooth
    gravity = '9.8 0 0'
    outputs = exodus
    output_properties = 'profile_mixture_density gas_velocity liquid_velocity void_fraction flow_pattern current_phase gas_density liquid_density mass_fraction density specific_heat temperature'
  [../]
[]

# [BCs]
#   [./pbc]
#     type = DirichletBC
#     variable = p
#     boundary = left
#     value = '2.07e6'
#   [../]
#   [./qbc]
#     type = DirichletBC
#     variable = q
#     boundary = right
#     value = 0.046296296
#   [../]
# []

[Variables]
  [./h]
    initial_condition =1500000
  [../]
  [./p]
    initial_condition = 2.07e6
  [../]
  [./q]
    scaling = 1e-2
    initial_condition = 0.046296296
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
  active = p2
  [./p1]
    type = SMP
    full = true
    #petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-pc_type -pc_hypre_type -snes_type -snes_linesearch_type -sub_pc_factor_shift_type'
    petsc_options_value = 'hypre boomeramg newtonls basic NONZERO'
  [../]
  [./p2]
    type = SMP
    full = true
    #petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-pc_type -sub_pc_type -snes_type -snes_linesearch_type -sub_pc_factor_shift_type -ksp_gmres_restart'
    petsc_options_value = 'asm lu newtonls basic NONZERO 51'
  [../]
  [./p3]
    type = SMP
    full = true
    #petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-pc_type -ksp_type -sub_pc_type -snes_type -snes_linesearch_type -pc_asm_overlap -sub_pc_factor_shift_type -ksp_gmres_restart'
    petsc_options_value = 'asm gmres lu newtonls basic 2 NONZERO 51'
  [../]
  [./p4]
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
