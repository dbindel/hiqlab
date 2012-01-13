-- Start clock
ts   = os.clock()

-- Loading the Lua input file
f0 = 0
runfile = loadfile('testg.lua')

-- Run the Lua input file
runfile()

-- Initialize the mesh
mpicomm = Get_Epetra_Comm()
mesh:initialize()
mesh:apply_bc()
numid  = mesh:get_numid()
mypid  = mpicomm:MyPID()
numproc= mpicomm:NumProc()
if (mypid==0) then
  print('Number of processes:',numproc)
  print("NDOF:",numid)
  -- Open file to write data
  filename = string.format("trilinos_results_%d.txt",numproc)
  f = io.open(filename,"w")
  f:write('Number of processes:',numproc,'\n')
  f:write("NDOF:",numid,'\n')
end

-- Assemble distributed stiffness
K       = mesh:assemble_dR_trilinos(1,0,0,1)

-- Assemble RHS
F = mesh:assemble_R_trilinos()
F:Scale(-1)

-- Assemble LHS
U = Epetra_Vector:new(numid)

-- Set up the linear system 
ELP= Epetra_LinearProblem:new(qCrs2RowMatrix(K),U,F)

-- Construct AztecOO
AZS = AztecOO:new(ELP)
AZS:SetAztecOption(AZ_solver, AZ_gmres_condnum)
AZS:SetAztecOption(AZ_output, 10)
AZS:SetAztecOption(AZ_kspace, 400)

-- Construc ML preconditioner
mlpl = ParameterList:new()
ML_Epetra_SetDefaults(mlpl, "SA")
mlpl:set_int(   "output", 0)
mlpl:set_string("prec type", "full-MGV")
mlpl:set_int(   "print unused",  1)
mlpl:set_int(   "max levels", 10)
mlpl:set_string("increasing or decreasing", "increasing")
mlpl:set_string("aggregation: type", "ParMETIS")
mlpl:set_string("smoother: type", "Aztec")
mlpl:set_double("smoother: damping factor", 0.67)
mlpl:set_int(   "smoother: sweeps", 2)
mlpl:set_string("smoother: pre or post", "both")
mlpl:set_string("coarse: type", "Amesos-KLU")

-- Construct Aztec Options and Params for ML smoother
Aop = AztecOptionsParams:new()
Aop:SetOption(AZ_precond, AZ_dom_decomp)
Aop:SetOption(AZ_subdomain_solve, AZ_ilu)
Aop:SetPL(mlpl)

MLPrec = qCreate_ML_Epetra_MultiLevelPreconditioner(K,mlpl)
AZS:SetPrecOperator(MLPrec)

-- Solve problem 
ts_s   = os.clock()
AZS:Iterate(2000, 1e-9)
te_s    = os.clock()
its     = AZS:NumIters()
resid_tr= AZS:TrueResidual()
resid_sr= AZS:ScaledResidual()
resid_rr= AZS:RecursiveResidual()

-- Check residual
R = Epetra_Vector:new(numid)
K:Multiply(false, U, R)
R:Update( -1, F, 1);

resid = {}
residf= {}
R:Norm2(resid)
F:Norm2(residf)
if (mypid==0) then
  print('Norm of residual(recu)  :',resid_rr)
  print('Norm of residual(abs)   :',resid[1])
  print('Norm of residual(rel)   :',resid[1]/residf[1])
  print('Number of Iterations    :',its)
  f:write('Norm of residual(recu)  :',resid_rr,'\n')
  f:write('Norm of residual(abs)   :',resid[1],'\n')
  f:write('Norm of residual(rel)   :',resid[1]/residf[1],'\n')
  f:write('Number of Iterations    :',its,'\n')
end

-- Print results
--[[
ToMatrixMarketFile('K.mm', K)
ToMatrixMarketFile('U.mm', U)
ToMatrixMarketFile('F.mm', F)
--]]

-- Write file for plotting
--[[
print('MyPID:',mypid)
Mesh_SetU_Epetra_MultiVector(mesh,U,0)
if mypid==0 then
  dxf = DXFile:new('testg0')
  dxf:writemesh(mesh)
  dxf:delete()
end
--]]

-- Delete objects
R:delete()
MLPrec:delete()
Aop:delete()
mlpl:delete()
AZS:delete()
ELP:delete()
U:delete()
F:delete()
K:delete()
mesh:delete()

-- Stop clock
te   = os.clock()
if (mypid==0) then
   print('Elapsed time            :',te-ts)
   print('Elapsed time(Only solve):',te_s-ts_s)
   f:write('Elapsed time            :',te-ts,'\n')
   f:write('Elapsed time(Only solve):',te_s-ts_s,'\n')
   f:close()
end
