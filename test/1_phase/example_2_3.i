# "Applied multiphase flow in pipes and flow assurance oil and gas production"
# Al-Safran, E., Brill, J. P., 2017
# Example 2.3: Determine the bottom hole pressure in the compressible gas well?

[Mesh]
  type = GeneratedMesh
  dim = 1
  xmin = 0
  xmax = 3048
  nx = 762
[]

[UserObjects]
  [./eos]
    type = MoskitoEOSIdealGas
    molar_mass = 21.723323e-3
    cp = 1000
  [../]
  [./viscosity]
    type = MoskitoViscosityConst
    viscosity = 0.000018
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
    well_diameter = 0.062
    roughness_type = rough
    roughness = 2.13e-5
    gravity = '9.8 0 0'
  [../]
[]

[BCs]
  [./pbc]
    type = DirichletBC
    variable = p
    boundary = left
    value = '13.83e6'
  [../]
  [./qbc]
    type = DirichletBC
    variable = q
    boundary = left
    value = 0
  [../]
[]

[Variables]
  [./h]
    [./InitialCondition]
      type = FunctionIC
      function = '(43.3+0.024606*x)*1000'
      variable = h
    [../]
  [../]
  [./p]
    initial_condition = 13.83e6
  [../]
  [./q]
    scaling = 1e-2
    initial_condition = 0
  [../]
[]

[Kernels]
  [./hkernel]
    type = NullKernel
    variable = h
  [../]
  [./pkernel]
    type = MoskitoMass1P
    variable = p
    flowrate = q
    enthalpy = h
  [../]
  [./qkernel]
    type = MoskitoMomentum1P
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
