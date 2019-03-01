[Mesh]
  type = FileMesh
  file = well_reservoir.msh
[]

# [MeshModifiers]
#   # [./b1]
#   #   type = BlockDeleter
#   #   block_id = 1
#   # [../]
#   [./b2]
#     type = BlockDeleter
#     block_id = 3
#   [../]
# []

[UserObjects]
  [./eos]
    type = MoskitoEOSIdealFluid
    bulk_modulus = 1e+9
  [../]
  [./viscosity]
    type = MoskitoViscosityConst
    viscosity = 0.001
  [../]
  [./rock_uo]
    type =  TigerPermeabilityConst
    permeability_type = isotropic
    k0 = '1.0e-14'
  [../]
[]

[Modules]
  [./FluidProperties]
    [./water_uo]
      type = TigerWaterConst
    [../]
  [../]
[]

[Materials]
  [./wells]
    type = MoskitoFluidWell1P
    pressure = p
    enthalpy = h
    flowrate = q
    well_direction = -z
    eos_uo = eos
    viscosity_uo = viscosity
    well_diameter = 0.152
    roughness_type = smooth
    gravity = '0 0 -9.8'
    block = 'injwell prowell'
    output_properties = 'density well_velocity'
    outputs = exodus
  [../]
  [./matrix_g]
    type = TigerGeometryMaterial
    porosity = 0.4
    block = matrix
  [../]
  [./matrix_f]
    type = TigerFluidMaterial
    fp_uo = water_uo
    block = matrix
  [../]
  [./matrix_h]
    type = TigerHydraulicMaterialH
    pressure = p
    compressibility = 1.0e-9
    kf_uo = rock_uo
    block = matrix
    output_properties = 'darcy_velocity'
    outputs = exodus
  [../]
[]

[BCs]
  # [./pbc_inj]
  #   type = DirichletBC
  #   variable = p
  #   boundary = injpoint
  #   value = 0
  # [../]
  # [./pbc_pro]
  #   type = DirichletBC
  #   variable = p
  #   boundary = propoint
  #   value = -0
  # [../]
  [./pbc_matrix]
    type = DirichletBC
    variable = p
    boundary = surround
    value = 15e6
  [../]
  [./qbc_inj]
    type = DirichletBC
    variable = q
    boundary = injpoint
    value = 0.03
  [../]
  [./qbc_pro]
    type = DirichletBC
    variable = q
    boundary = propoint
    value = -0.03
  [../]
[]

[Variables]
  [./h]
    initial_condition = 0
    block = 'injwell prowell'
  [../]
  [./p]
    [./InitialCondition]
      type = FunctionIC
      variable = p
      function = -1000*9.8*z
    [../]
    scaling = 1e-6
  [../]
  [./q]
    initial_condition = 0
    block = 'injwell prowell'
  [../]
[]

[Kernels]
  [./hkernel]
    type = NullKernel
    variable = h
    block = 'injwell prowell'
  [../]
  [./pkernelm]
    type = TigerHydraulicKernelH
    variable = p
    block = matrix
  [../]
  [./pkernelw]
    type = MoskitoMass
    variable = p
    flowrate = q
    enthalpy = h
    block = 'injwell prowell'
  [../]
  [./qkernel]
    type = MoskitoMomentum
    variable = q
    pressure = p
    enthalpy = h
    block = 'injwell prowell'
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
