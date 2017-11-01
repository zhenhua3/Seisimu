# function test_damping()

  include("Project/MODL_UNSPLT_PML/src/ModelSetup/InitWavefield2D.jl")
  include("Project/MODL_UNSPLT_PML/src/ModelSetup/Medium.jl")
  include("Project/MODL_UNSPLT_PML/src/ModelSetup/PMLDamping.jl")

  v = Float64[3000, 4000, 1700, 2200]
  rho = [2.5, 2.5]
  LN = 4
  Depth = 2000.0
  Horiz = 3000.0
  PF = 120.0
  ext = 10
  iflag = 2

  medium = ModInit()
  medium = model(v, rho, LN, Depth, Horiz, PF, medium, ext, iflag)

  WF = InitWF2D(medium.NDep, medium.NHor, ext, iflag)

  dampbound = DampBound(WF, medium, ext, iflag)
# end

# using PyCall
# @pyimport matplotlib.pyplot as plt
# plt.imshow(tmpdVxx)
# plt.show()
