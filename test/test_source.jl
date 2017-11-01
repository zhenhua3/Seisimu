# function test_Source()
  include("Project/MODL_UNSPLT_PML/src/ModelSetup/Medium.jl")
  include("Project/MODL_UNSPLT_PML/src/ModelSetup/InitWavefield2D.jl")
  include("Project/MODL_UNSPLT_PML/src/ModelSetup/Source/Source.jl")
  include("Project/MODL_UNSPLT_PML/src/ModelSetup/Source/Addsource.jl")
  include("Project/MODL_UNSPLT_PML/src/ModelSetup/Source/Wavelet.jl")

  v = Float64[3000, 4000, 1700, 2200]
  rho = [2.5, 2.5]
  LN = 4
  Depth = 2000
  Horiz = 3000.0
  PF = [20,30,40]
  maxPF = maximum(PF)
  ext = 10
  iflag = 1

  medium = ModInit()
  medium = model(v, rho, LN, Depth, Horiz, maxPF, medium, ext, iflag)
  WF = InitWF2D(medium.NDep, medium.NHor, ext, iflag)

  # test single source
  position = [300.0;400]
  sot = 0.03
  stp = "expl"

  waveform = Ricker(PF[1]/2,medium.dt)

  T = 1.0

  t_source = source(position,sot,ext,stp,medium,waveform,T,WF,iflag)

  # test multiple source
  sn = 3
  physour = Array{MultiSource}(sn)
  waveform1 = Ricker(20/2,medium.dt)
  physour[1] = MultiSource("sfx",[300;400],20,0.1,waveform1)

  waveform2 = Ricker(30/2,medium.dt)
  physour[2] = MultiSource("sfz",[400;400],30,0.1,waveform2)

  waveform3 = Ricker(40/2,medium.dt)
  physour[3] = MultiSource("DC",[100;400],40,0.1,waveform3)

  source(sn, physour, ext, medium, T, WF, iflag)
  # return t_source
# end
