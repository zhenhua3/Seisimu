function test_Addsource()

  include("Project/MODL_UNSPLT_PML/src/ModelSetup/Medium.jl")
  include("Project/MODL_UNSPLT_PML/src/ModelSetup/InitWavefield2D.jl")
  include("Project/MODL_UNSPLT_PML/src/ModelSetup/Source/Source.jl")
  include("Project/MODL_UNSPLT_PML/src/ModelSetup/Source/Addsource.jl")
  include("Project/MODL_UNSPLT_PML/src/ModelSetup/Source/Wavelet.jl")
  include("Project/MODL_UNSPLT_PML/src/ModelSetup/Source/Addsource.jl")

  v = Float64[3000, 4000, 1700, 2200]
  rho = [2.5, 2.5]
  LN = 4
  Depth = 2000.0
  Horiz = 3000.0
  PF = [120.0,80, 60]
  maxPF = maximum(PF)
  ext = 10
  iflag = 1

  medium = ModInit()
  medium = model(v, rho, LN, Depth,Horiz, maxPF, medium, ext, iflag)

  WF = InitWF2D(medium.PML_NDep, medium.PML_NHor, ext, iflag)

  sn = 3
  position = ([300.0,400], [500.0, 700], [600.0,800])
  ot = [0.03, 0.1, 0.05]
  tp = ["expl", "sfx", "sfz"]


  waveform1 = Ricker(PF[1]/2,medium.dt)
  waveform2 = Ricker(PF[2]/2,medium.dt)
  waveform3 = Ricker(PF[3]/2,medium.dt)
  waveform = (waveform1, waveform2, waveform3)

  T = 3.0

  t_source = source(sn,position,ot,ext,tp,medium,waveform,T, iflag)

for i = 1:1000
  Addsource!(WF, sn, t_source, i)

println(WF.VecBDTxx1[t_source[1].BDposn[1,1]])
end

end
