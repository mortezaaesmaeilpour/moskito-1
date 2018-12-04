[Mesh]
  type = GeneratedMesh
  dim = 1
  xmin = 0
  xmax = 1000
  nx = 1000
[]

[MeshModifiers]
  [./block1]
    type = SubdomainBoundingBox
    block_id = 1
    bottom_left = '-1 -1 -1'
    top_right = '500 1 1'
  [../]
[]

[UserObjects]
  [./eos]
    type = MoskitoEOSIdealFluid
    reference_pressure = 0
    bulk_modulus = 2e+09
  [../]
  [./viscosity]
    type = MoskitoViscosityConst
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
    well_diameter = 0.2
    roughness_type = smooth
    output_properties = 'density'
    gravity = '10 0 0'
    block = 0
  [../]
  [./area1]
    type = MoskitoFluidWell1P
    pressure = p
    enthalpy = h
    flowrate = q
    well_direction = x
    eos_uo = eos
    viscosity_uo = viscosity
    well_diameter = 0.25
    roughness_type = smooth
    output_properties = 'density'
    gravity = '10 0 0'
    block = 1
  [../]
[]

[BCs]
  [./pbc]
    type = DirichletBC
    variable = p
    boundary = left
    value = 0
  [../]
  [./qbc]
    type = DirichletBC
    variable = q
    boundary = right
    value = 0
  [../]
[]

[Variables]
  [./h]
    initial_condition = 83950
  [../]
  [./p]
    initial_condition = 0
  [../]
  [./q]
    scaling = 1e-3
    initial_condition = 0
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

[Outputs]
  [./out]
    type = Exodus
    output_material_properties = true
  [../]
  [./out1]
    type = VariableResidualNormsDebugOutput
  [../]
[]
