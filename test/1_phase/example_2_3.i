# "Applied multiphase flow in pipes and flow assurance oil and gas production"
# Al-Safran, E., Brill, J. P., 2017
# Example 2.3: Determine the bottom hole pressure in the compressible gas well?

[Mesh]
  type = GeneratedMesh
  dim = 1
  xmin = 0
  xmax = 3048
  nx = 1524
  allow_renumbering = false
[]

[UserObjects]
  [./eos]
    type = MoskitoEOSNaturalGas
    # type = MoskitoEOSIdealGas
    molar_mass = 2.16e-2
    specific_gravity = 0.75
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
    gravity = '9.81 0 0'
    # manual_friction_factor = 0.0149
    output_properties = 'density well_velocity well_reynolds_no well_moody_friction temperature drho_dp drho_dT hx px'
  [../]
[]

[AuxVariables]
  [./rho]
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./rho]
    type = MaterialRealAux
    property = density
    variable = rho
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
    value = 0.0093034
  [../]
[]

[Variables]
  [./h]
    [./InitialCondition]
      type = FunctionIC
      function = '(273.15+43.3+0.024606*x)*1000'
      variable = h
    [../]
  [../]
  [./p]
    initial_condition = 13.83e6
  [../]
  [./q]
    scaling = 1e-4
    initial_condition = 0.0093034
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
  nl_rel_tol = 1e-9
  nl_abs_tol = 1e-9
  nl_max_its = 50
  solve_type = NEWTON
[]

[Postprocessors]
  [./p0]
    type = NodalVariableValue
    variable = p
    nodeid = 1
  [../]
  [./p1]
    type = NodalVariableValue
    variable = p
    nodeid = 381
  [../]
  [./p2]
    type = NodalVariableValue
    variable = p
    nodeid = 1143
  [../]
  [./p3]
    type = NodalVariableValue
    variable = p
    nodeid = 1524
  [../]
  [./q0]
    type = NodalVariableValue
    variable = q
    nodeid = 1
    scale_factor = 331.227769182
  [../]
  [./q1]
    type = NodalVariableValue
    variable = q
    nodeid = 381
    scale_factor = 331.227769182
  [../]
  [./q2]
    type = NodalVariableValue
    variable = q
    nodeid = 1143
    scale_factor = 331.227769182
  [../]
  [./q3]
    type = NodalVariableValue
    variable = q
    nodeid = 1524
    scale_factor = 331.227769182
  [../]
  [./rho0]
    type = ElementalVariableValue
    variable = rho
    elementid = 1
  [../]
  [./rho1]
    type = ElementalVariableValue
    variable = rho
    elementid = 381
  [../]
  [./rho2]
    type = ElementalVariableValue
    variable = rho
    elementid = 1143
  [../]
[]

[Outputs]
  # exodus = true
  [./out]
    type = Exodus
    output_material_properties = true
  [../]
  # [./res]
  #   type = VariableResidualNormsDebugOutput
  #   nonlinear_residuals = false
  # [../]
[]
