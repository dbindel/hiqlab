-- Loading the Lua input file
f0 = 0
order = 1;
dense = 5;
runfile = loadfile('plate3d.lua')
--runfile = loadfile('pml1d.lua')

-- Run the Lua input file
runfile()

-- Initialize the mesh
mesh:initialize()
mesh:apply_bc()
numid = mesh:get_numid()

-- Assemble distributed stiffness
K = mesh:assemble_dR_trilinos(1,0,0,1)
M = mesh:assemble_dR_trilinos(0,0,1,1)

-- Set parameters
PL = ParameterList:new()
PL:set_int   ("Solver", Anasazi_BlockKrylovSchur)
PL:set_int   ("Block Size",        2)
PL:set_int   ("Max Blocks",        4)
PL:set_int   ("Max Restarts",    100)
PL:set_double("Convergence Tolerance",          1.0e-9)
PL:set_string("Sigma",          "LM")    -- Sort manager parameters
PL:set_int("Verbosity", Anasazi_Warning+Anasazi_FinalSummary+Anasazi_IterationDetails+Anasazi_Debug+Anasazi_OrthoDetails) -- Output manager parameters
--PL:set_int("Verbosity", Anasazi_Warning) -- Output manager parameters

-- Set memory for eigenvector
PL:set_int("IsSymmetric", 0)
nev= 2
Vr = Epetra_MultiVector:new(numid,1*nev)

-- Compute eigenvalues
dr = {}
di = {}
w0 = 0
status = compute_eigs_anasazi1(
K, 
M, 
nev, 
PL, dr, di, Vr)
undo_spectral_trans(nev,w0,0,1,dr,di)

-- Print eigenvalues
for i = 1,nev do
   print(dr[i],' ',di[i])
end

-- Print results
--[[
ToMatrixMarketFile('K.mm', K)
ToMatrixMarketFile('M.mm', M)
ToMatrixMarketFile('Vr.mm',Vr)
ToMatrixMarketFile('Vi.mm',Vi)
--]]

-- Delete objects
Vr:delete()
PL:delete()
M:delete()
K:delete()
mesh:delete()
