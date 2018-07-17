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
  [../]
[]

[Materials]
  [./area0]
    type = MoskitoCMomentumMaterial
    density = rho
    eos_UO = eos
    well_diameter = 0.2
    block = 0
  [../]
  [./area1]
    type = MoskitoCMomentumMaterial
    well_diameter = 0.25
    density = rho
    eos_UO = eos
    block = 1
  [../]
[]

[AuxVariables]
  [./p]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./area]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./pkernel]
    type = MaterialRealAux
    variable = p
    property = 'pressure'
  [../]
  [./areakernel]
    type = MaterialRealAux
    variable = area
    property = 'well_area'
  [../]
[]

[BCs]
  [./rhobc]
    type = DirichletBC
    variable = rho
    boundary = left
    value = 1000
  [../]
  [./qbc]
    type = DirichletBC
    variable = q
    boundary = left
    value = 0.05
  [../]
[]

[Variables]
  [./rho]
    initial_condition = 1000
  [../]
  [./q]
    scaling = 1e-6
    initial_condition = 0.05
  [../]
[]

[Kernels]
  [./rhokernel]
    type = MoskitoCMass
    variable = rho
    vol_flow_rate = q
  [../]
  [./qkernel]
    type = MoskitoCMomentum
    variable = q
    density = rho
  [../]
[]

[Preconditioning]
  active = 'p3'
  [./p1]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_hypre_type'
    petsc_options_value = 'hypre boomeramg'
  [../]
  [./p2]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type -ksp_gmres_restart'
    petsc_options_value = 'asm lu NONZERO 51'
  [../]
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
  nl_abs_tol = 1e-10
  nl_max_its = 50
  solve_type = NEWTON
[]

[Postprocessors]
  [./rho]
    type = VariableResidual
    variable = rho
  [../]
  [./q]
    type = VariableResidual
    variable = q
  [../]
[]

[Outputs]
  exodus = true
  print_linear_residuals = true
  # [./test]
  #   type = Exodus
  #   output_material_properties = true
  # [../]
[]
