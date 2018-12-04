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
    type = MoskitoEOSIdealFluid
    bulk_modulus = 2e+012
  [../]
  [./viscosity]
    type = MoskitoViscosityConst
    viscosity = 0.001
  [../]
[]

[Materials]
  [./area0]
    type = MoskitoFluidWell1P
    pressure = p
    enthalpy = h
    flowrate = q
    well_direction = x
    eos_uo = eos
    viscosity_uo = viscosity
    well_diameter = 0.152
    roughness_type = smooth
    gravity = '0 -9.8 0'
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
    value = -0.046296296
  [../]
[]

[Variables]
  [./h]
    initial_condition = 0
  [../]
  [./p]
    initial_condition = 2.07e6
  [../]
  [./q]
    scaling = 1e-2
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
