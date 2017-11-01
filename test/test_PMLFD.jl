# function test_PMLFD()

  include("Project/MODL_UNSPLT_PML/src/ModelSetup/InitWavefield2D.jl")
  include("Project/MODL_UNSPLT_PML/src/ModelSetup/Medium.jl")
  include("Project/MODL_UNSPLT_PML/src/ModelSetup/AddPML/PMLDamping.jl")
  include("Project/MODL_UNSPLT_PML/src/FiniteDifference/2D/FDCoeff.jl")
  include("Project/MODL_UNSPLT_PML/src/FiniteDifference/2D/Dx2D.jl")
  include("Project/MODL_UNSPLT_PML/src/FiniteDifference/2D/Dz2D.jl")
  include("Project/MODL_UNSPLT_PML/src/FiniteDifference/2D/FDMtx.jl")
  include("Project/MODL_UNSPLT_PML/src/FiniteDifference/2D/Onestep2D.jl")
  include("Project/MODL_UNSPLT_PML/src/ModelSetup/AddPML/PMLFD.jl")

  v = Float64[3000, 4000, 1700, 2200]
  rho = [2.5, 2.5]
  LN = 4
  Depth = 2000.0
  Horiz = 3000.0
  PF = 20.0
  ext = 10
  iflag = 2

  medium = ModInit()
  medium = model(v, rho, LN, Depth, Horiz, PF, medium, ext, iflag)

  WF = InitWF2D(medium.NDep, medium.NHor, ext, iflag)

  dampbound = DampBound(WF, medium, ext, iflag)

  FDC = FDCoeff(4)

  FD = Fdmtx(WF, medium, FDC, ext)

  pml = AddPML(dampbound, FD)

  AddPML!(pml, dampbound, FD)

  PMLFD!(pml, FD)

# end
