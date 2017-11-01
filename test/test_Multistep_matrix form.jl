using PyCall
@pyimport matplotlib.pyplot as plt

include("Project/MODL_UNSPLT_PML/src/ModelSetup/InitWavefield2D.jl")
include("Project/MODL_UNSPLT_PML/src/ModelSetup/Medium.jl")
include("Project/MODL_UNSPLT_PML/src/FiniteDifference/2D/FDCoeff.jl")
include("Project/MODL_UNSPLT_PML/src/FiniteDifference/2D/Dx2D.jl")
include("Project/MODL_UNSPLT_PML/src/FiniteDifference/2D/Dz2D.jl")
include("Project/MODL_UNSPLT_PML/src/FiniteDifference/2D/FDMtx.jl")
include("Project/MODL_UNSPLT_PML/src/ModelSetup/AddPML/PMLDamping.jl")
include("Project/MODL_UNSPLT_PML/src/ModelSetup/AddPML/PMLFD.jl")
include("Project/MODL_UNSPLT_PML/src/FiniteDifference/2D/Onestep2D.jl")
include("Project/MODL_UNSPLT_PML/src/ModelSetup/Source/Source.jl")
include("Project/MODL_UNSPLT_PML/src/ModelSetup/Source/Addsource.jl")
include("Project/MODL_UNSPLT_PML/src/ModelSetup/Source/Wavelet.jl")
include("Project/MODL_UNSPLT_PML/src/FiniteDifference/2D/MultiStep2D_nPML.jl")
include("Project/MODL_UNSPLT_PML/src/DataIO/WriteData.jl")
include("Project/MODL_UNSPLT_PML/src/DataIO/WriteInfo.jl")


  ############## ModelSetup ###############
  vp = 3000.0
  vs = 1700.0
  rho = 2500.0
  Depth = 1000.0
  Horiz = 1000.0
  PF = 80.0
  ext = 20

  iflag = 1 # free surface

  medium = ModInit()
  medium = model(vp, vs, rho, Depth, Horiz, PF, medium, ext, iflag)

  ############### WF initialization #############
  WF = InitWF2D(medium.NDep, medium.NHor, ext, iflag)

  ############### Finite Difference Matrix ################
  FDC = FDCoeff(4)

  FD = Fdmtx(WF, medium, FDC, ext)

  ############ Source ############
  sn = 1
  position = [500.0,500]
  ot = 0.005
  tp = "expl"

  waveform1 = Ricker(PF[1]/2,medium.dt)
  waveform = (waveform1)

  T = 0.15

  t_source = source(sn,position,ot,ext,tp,medium,waveform,T, iflag)


  SpatialSize = 822057
  ### reorgnize finite differential matrix
  L = spzeros(SpatialSize,SpatialSize)
  (nz_VxBTxx,nx_VxBTxx) = size(FD.VxBTxx)
  (nz_VxBTxz,nx_VxBTxz) = size(FD.VxBTxz)
  (nz_VzBTzz,nx_VzBTzz) = size(FD.VzBTzz)
  (nz_VzBTxz,nx_VzBTxz) = size(FD.VzBTxz)
  (nz_TxxBVx,nx_TxxBVx) = size(FD.TxxBVx)
  (nz_TxxBVz,nx_TxxBVz) = size(FD.TxxBVz)
  (nz_TzzBVx,nx_TzzBVx) = size(FD.TzzBVx)
  (nz_TzzBVz,nx_TzzBVz) = size(FD.TzzBVz)
  (nz_TxzBVx,nx_TxzBVx) = size(FD.TxzBVx)
  (nz_TxzBVz,nx_TxzBVz) = size(FD.TxzBVz)

  L[1:nz_VxBTxx,(nx_TxxBVx+nx_TxxBVz+1):(nx_TxxBVx+nx_TxxBVz+nx_VxBTxx)] = FD.VxBTxx
  L[1:nz_VxBTxz,(nx_TxxBVx+nx_TxxBVz+nx_VxBTxx+nx_VzBTzz+1):(nx_TxxBVx+nx_TxxBVz+nx_VxBTxx+nx_VzBTzz+nx_VxBTxz)] = FD.VxBTxz

  L[nz_VxBTxx+1:nz_VxBTxx+nz_VzBTzz,(nx_TxxBVx+nx_TxxBVz+nx_VxBTxx+1):(nx_TxxBVx+nx_TxxBVz+nx_VxBTxx+nx_VzBTzz)] = FD.VzBTzz
  L[nz_VxBTxx+1:nz_VxBTxx+nz_VzBTzz,(nx_TxxBVx+nx_TxxBVz+nx_VxBTxx+nx_VzBTzz+1):(nx_TxxBVx+nx_TxxBVz+nx_VxBTxx+nx_VzBTzz+nx_VxBTxz)] = FD.VzBTxz

  L[nz_VxBTxx+nz_VzBTzz+1:nz_VxBTxx+nz_VzBTzz+nz_TxxBVx,1:nx_TxxBVx] = FD.TxxBVx
  L[nz_VxBTxx+nz_VzBTzz+1:nz_VxBTxx+nz_VzBTzz+nz_TxxBVx,nx_TxxBVx+1:nx_TxxBVx+nx_TxxBVz] = FD.TxxBVz
  L[nz_VxBTxx+nz_VzBTzz+1:nz_VxBTxx+nz_VzBTzz+nz_TxxBVx,nx_TxxBVx+nx_TxxBVz+1:nx_TxxBVx+nx_TxxBVz+nx_VxBTxx] = FD.TxxBVx*FD.VxBTxx
  L[nz_VxBTxx+nz_VzBTzz+1:nz_VxBTxx+nz_VzBTzz+nz_TxxBVx,nx_TxxBVx+nx_TxxBVz+nx_VxBTxx+1:nx_TxxBVx+nx_TxxBVz+nx_VxBTxx+nx_VzBTzz] = FD.TxxBVz*FD.VzBTzz
  L[nz_VxBTxx+nz_VzBTzz+1:nz_VxBTxx+nz_VzBTzz+nz_TxxBVx,nx_TxxBVx+nx_TxxBVz+nx_VxBTxx+nx_VzBTzz+1:nx_TxxBVx+nx_TxxBVz+nx_VxBTxx+nx_VzBTzz+nx_VxBTxz] = FD.TxxBVx*FD.VxBTxz+FD.TxxBVz*FD.VzBTxz

  L[nz_VxBTxx+nz_VzBTzz+nz_TxxBVx+1:nz_VxBTxx+nz_VzBTzz+nz_TxxBVx+nz_TzzBVx,1:nx_TxxBVx] = FD.TzzBVx
  L[nz_VxBTxx+nz_VzBTzz+nz_TxxBVx+1:nz_VxBTxx+nz_VzBTzz+nz_TxxBVx+nz_TzzBVx,nx_TxxBVx+1:nx_TxxBVx+nx_TxxBVz] = FD.TzzBVz
  L[nz_VxBTxx+nz_VzBTzz+nz_TxxBVx+1:nz_VxBTxx+nz_VzBTzz+nz_TxxBVx+nz_TzzBVx,nx_TxxBVx+nx_TxxBVz+1:nx_TxxBVx+nx_TxxBVz+nx_VxBTxx] = FD.TzzBVx*FD.VxBTxx
  L[nz_VxBTxx+nz_VzBTzz+nz_TxxBVx+1:nz_VxBTxx+nz_VzBTzz+nz_TxxBVx+nz_TzzBVx,nx_TxxBVx+nx_TxxBVz+nx_VxBTxx+1:nx_TxxBVx+nx_TxxBVz+nx_VxBTxx+nx_VzBTzz] = FD.TzzBVz*FD.VzBTzz
  L[nz_VxBTxx+nz_VzBTzz+nz_TxxBVx+1:nz_VxBTxx+nz_VzBTzz+nz_TxxBVx+nz_TzzBVx,nx_TxxBVx+nx_TxxBVz+nx_VxBTxx+nx_VzBTzz+1:nx_TxxBVx+nx_TxxBVz+nx_VxBTxx+nx_VzBTzz+nx_VxBTxz] = FD.TzzBVx*FD.VxBTxz+FD.TzzBVz*FD.VzBTxz

  L[nz_VxBTxx+nz_VzBTzz+nz_TxxBVx+nz_TzzBVx+1:nz_VxBTxx+nz_VzBTzz+nz_TxxBVx+nz_TzzBVx+nz_TxzBVx,1:nx_TxxBVx] = FD.TxzBVx
  L[nz_VxBTxx+nz_VzBTzz+nz_TxxBVx+nz_TzzBVx+1:nz_VxBTxx+nz_VzBTzz+nz_TxxBVx+nz_TzzBVx+nz_TxzBVx,nx_TxxBVx+1:nx_TxxBVx+nx_TxxBVz] = FD.TxzBVz
  L[nz_VxBTxx+nz_VzBTzz+nz_TxxBVx+nz_TzzBVx+1:nz_VxBTxx+nz_VzBTzz+nz_TxxBVx+nz_TzzBVx+nz_TxzBVx,nx_TxxBVx+nx_TxxBVz+1:nx_TxxBVx+nx_TxxBVz+nx_VxBTxx] = FD.TxzBVx*FD.VxBTxx
  L[nz_VxBTxx+nz_VzBTzz+nz_TxxBVx+nz_TzzBVx+1:nz_VxBTxx+nz_VzBTzz+nz_TxxBVx+nz_TzzBVx+nz_TxzBVx,nx_TxxBVx+nx_TxxBVz+nx_VxBTxx+1:nx_TxxBVx+nx_TxxBVz+nx_VxBTxx+nx_VzBTzz] = FD.TxzBVz*FD.VzBTzz
  L[nz_VxBTxx+nz_VzBTzz+nz_TxxBVx+nz_TzzBVx+1:nz_VxBTxx+nz_VzBTzz+nz_TxxBVx+nz_TzzBVx+nz_TxzBVx,nx_TxxBVx+nx_TxxBVz+nx_VxBTxx+nx_VzBTzz+1:nx_TxxBVx+nx_TxxBVz+nx_VxBTxx+nx_VzBTzz+nx_VxBTxz] = FD.TxzBVx*FD.VxBTxz+FD.TxzBVz*FD.VzBTxz

  t1 = [WF.VecBDVx1;WF.VecBDVz1;WF.VecBDTxx1;WF.VecBDTzz1;WF.VecBDTxz1]

  t3 = [WF.VecBDVx3;WF.VecBDVz3;WF.VecBDTxx3;WF.VecBDTzz3;WF.VecBDTxz3]

  for it = 1:t_source.Tn

    Addsource!(WF, t_source, it)
    t1 = [WF.VecBDVx1;WF.VecBDVz1;WF.VecBDTxx1;WF.VecBDTzz1;WF.VecBDTxz1]

    t3 = t1 + L*t1

    WF.VecBDVx1 = t3[1:396*415]
    WF.VecBDVz1 = t3[396*415+1:396*415+395*416]
    WF.VecBDTxx1 = t3[396*415+395*416+1:396*415+395*416+396*416]
    WF.VecBDTzz1 = t3[396*415+395*416+396*416+1:396*415+395*416+2*396*416]
    WF.VecBDTxz1 = t3[396*415+395*416+2*396*416+1:396*415+395*416+2*396*416+395*415]

    if mod(it,100) == 0
    println(it)
    image = reshape(full(WF.VecBDTxx1+WF.VecBDTzz1),WF.N_BDTpp[1], WF.N_BDTpp[2])
    # image = reshape(t3[1:396*415],396,415)
    plt.show(plt.imshow(image))
    end

  end
