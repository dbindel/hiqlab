
-- Utility functions on top of PETSc --


function petsc_vector_dump(x, fname)
  local viewer = PetscViewerASCIIOpen(fname)
  PetscViewerSetFormat(viewer, "matlab")
  VecView(x, viewer)
  PetscViewerFlush(viewer)
  PetscViewerDestroy(viewer)
end


function petsc_vector_bdump(x, fname)
  local viewer = PetscViewerBinaryOpen(fname, 'FILE_MODE_WRITE')
  VecView(x, viewer)
  PetscViewerFlush(viewer)
  PetscViewerDestroy(viewer)
end


function petsc_vector_load(fname)
  local viewer = PetscViewerBinaryOpen(fname, 'FILE_MODE_READ')
  local vec = VecLoad(viewer, 'mpi')
  PetscViewerDestroy(viewer)
  return vec
end


function petsc_matrix_dump(x, fname)
  local viewer = PetscViewerASCIIOpen(fname)
  PetscViewerSetFormat(viewer, "matlab")
  MatView(x, viewer)
  PetscViewerFlush(viewer)
  PetscViewerDestroy(viewer)
end


function petsc_matrix_bdump(x, fname)
  local viewer = PetscViewerBinaryOpen(fname, 'FILE_MODE_WRITE')
  MatView(x, viewer)
  PetscViewerFlush(viewer)
  PetscViewerDestroy(viewer)
end


function petsc_matrix_load(fname)
  local viewer = PetscViewerBinaryOpen(fname, 'FILE_MODE_READ')
  local mat = MatLoad(viewer, 'seqbaij')
  PetscViewerDestroy(viewer)
  return mat
end


function petsc_solve_x(p)

  local A    = p[1] or p.A or error('Must have valid matrix')
  local b    = p[2] or p.b or error('Must have valid rhs')
  local opt  = p[3] or p.opt or {}

  local B    = p.B or opt.B or A
  local mesh = p.mesh or opt.mesh
  local rtol = p.rtol or opt.rtol or PETSC_DEFAULT
  local atol = p.atol or opt.atol or PETSC_DEFAULT
  local dtol = p.dtol or opt.dtol or PETSC_DEFAULT
  local maxits = p.maxits or opt.maxits or PETSC_DEFAULT

  local mA, nA = MatGetSize(A)
  local mB, nB = MatGetSize(B)
  local mb     = VecGetSize(b)
  assert(mA == nA and mA == mB and mB == nB and mB == mb)

  local x = VecCreate()
  VecSetSizes(x, PETSC_DECIDE, nA)
  VecSetFromOptions(x)

  local ksp = KSPCreate()
  KSPSetOperators(ksp, A, B, "DIFFERENT_NONZERO_PATTERN")
  KSPSetTolerances(ksp, rtol, atol, dtol, maxits)
  KSPSetFromOptions(ksp)
  local pc = KSPGetPC(ksp)
  PCSetCoordinatesFromMesh(pc, mesh)

  KSPSolve(ksp,b,x)
  -- KSPView(ksp, PETSC_VIEWER_STDOUT_WORLD)
  local rnorm = KSPGetResidualNorm(ksp)
  local its   = KSPGetIterationNumber(ksp)
  KSPDestroy(ksp)

  return x, rnorm, its
end
