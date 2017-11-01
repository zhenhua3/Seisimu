function Onestep2D!(WF1::WF2D, WF2::WF2D, FD::NSFDMtx, pml::dampCoef)
# For unsplit PML boundary simulation

  VxBTxx = FD.VxBTxx*WF1.VecBDTxx
  VxBTxz = FD.VxBTxz*WF1.VecBDTxz
  VzBTzz = FD.VzBTzz*WF1.VecBDTzz
  VzBTxz = FD.VzBTxz*WF1.VecBDTxz

  pml.PVxBTxx = (pml.bxVx .* pml.PVxBTxx + pml.PVxBTxx + pml.axVx .* VxBTxx)
  pml.PVxBTxz = (pml.bzVx .* pml.PVxBTxz + pml.PVxBTxz + pml.azVx .* VxBTxz)
  pml.PVzBTzz = (pml.bzVz .* pml.PVzBTzz + pml.PVzBTzz + pml.azVz .* VzBTzz)
  pml.PVzBTxz = (pml.bxVz .* pml.PVzBTxz + pml.PVzBTxz + pml.axVz .* VzBTxz)

  WF2.VecBDVx = WF1.VecBDVx + VxBTxx + VxBTxz + pml.PVxBTxx + pml.PVxBTxz
  WF2.VecBDVz = WF1.VecBDVz + VzBTzz + VzBTxz + pml.PVzBTzz + pml.PVzBTxz

  TxxBVx = FD.TxxBVx*WF2.VecBDVx
  TxxBVz = FD.TxxBVz*WF2.VecBDVz
  TzzBVx = FD.TzzBVx*WF2.VecBDVx
  TzzBVz = FD.TzzBVz*WF2.VecBDVz
  TxzBVx = FD.TxzBVx*WF2.VecBDVx
  TxzBVz = FD.TxzBVz*WF2.VecBDVz

  pml.PTxxBVx = (pml.bxTxx .* pml.PTxxBVx + pml.PTxxBVx + pml.axTxx .* TxxBVx)
  pml.PTxxBVz = (pml.bzTxx .* pml.PTxxBVz + pml.PTxxBVz + pml.azTxx .* TxxBVz)
  pml.PTzzBVx = (pml.bxTzz .* pml.PTzzBVx + pml.PTzzBVx + pml.axTzz .* TzzBVx)
  pml.PTzzBVz = (pml.bzTzz .* pml.PTzzBVz + pml.PTzzBVz + pml.azTzz .* TzzBVz)
  pml.PTxzBVx = (pml.bzTxz .* pml.PTxzBVx + pml.PTxzBVx + pml.azTxz .* TxzBVx)
  pml.PTxzBVz = (pml.bxTxz .* pml.PTxzBVz + pml.PTxzBVz + pml.axTxz .* TxzBVz)

  WF2.VecBDTxx = WF1.VecBDTxx + TxxBVx + TxxBVz + pml.PTxxBVx + pml.PTxxBVz
  WF2.VecBDTzz = WF1.VecBDTzz + TzzBVx + TzzBVz + pml.PTzzBVx + pml.PTzzBVz
  WF2.VecBDTxz = WF1.VecBDTxz + TxzBVx + TxzBVz + pml.PTxzBVx + pml.PTxzBVz

  WF1.VecBDVx = copy(WF2.VecBDVx)
  WF1.VecBDVz = copy(WF2.VecBDVz)
  WF1.VecBDTxx = copy(WF2.VecBDTxx)
  WF1.VecBDTzz = copy(WF2.VecBDTzz)
  WF1.VecBDTxz = copy(WF2.VecBDTxz)

end

function Onestep2D!(WF1::WF2D, WF2::WF2D, FD::SFDMtx)
  WF2.VecBDVxx = FD.VxxBVxx*WF1.VecBDVxx + FD.VxxBTxx*(WF1.VecBDTxxx+WF1.VecBDTxxz)
  WF2.VecBDVxz = FD.VxzBVxz*WF1.VecBDVxz + FD.VxzBTxz*(WF1.VecBDTxzx+WF1.VecBDTxzz)
  WF2.VecBDVzx = FD.VzxBVzx*WF1.VecBDVzx + FD.VzxBTxz*(WF1.VecBDTxzx+WF1.VecBDTxzz)
  WF2.VecBDVzz = FD.VzzBVzz*WF1.VecBDVzz + FD.VzzBTzz*(WF1.VecBDTzzx+WF1.VecBDTzzz)

  WF2.VecBDTxxx = FD.TxxxBTxxx*WF1.VecBDTxxx + FD.TxxxBVx*(WF2.VecBDVxx+WF2.VecBDVxz)
  WF2.VecBDTxxz = FD.TxxzBTxxz*WF1.VecBDTxxz + FD.TxxzBVz*(WF2.VecBDVzx+WF2.VecBDVzz)
  WF2.VecBDTzzx = FD.TzzxBTzzx*WF1.VecBDTzzx + FD.TzzxBVx*(WF2.VecBDVxx+WF2.VecBDVxz)
  WF2.VecBDTzzz = FD.TzzzBTzzz*WF1.VecBDTzzz + FD.TzzzBVz*(WF2.VecBDVzx+WF2.VecBDVzz)
  WF2.VecBDTxzx = FD.TxzxBTxzx*WF1.VecBDTxzx + FD.TxzxBVz*(WF2.VecBDVzx+WF2.VecBDVzz)
  WF2.VecBDTxzz = FD.TxzzBTxzz*WF1.VecBDTxzz + FD.TxzzBVx*(WF2.VecBDVxx+WF2.VecBDVxz)

  WF1.VecBDVxx = copy(WF2.VecBDVxx)
  WF1.VecBDVxz = copy(WF2.VecBDVxz)
  WF1.VecBDVzx = copy(WF2.VecBDVzx)
  WF1.VecBDVzz = copy(WF2.VecBDVzz)
  WF1.VecBDTxxx = copy(WF2.VecBDTxxx)
  WF1.VecBDTxxz = copy(WF2.VecBDTxxz)
  WF1.VecBDTzzx = copy(WF2.VecBDTzzx)
  WF1.VecBDTzzz = copy(WF2.VecBDTzzz)
  WF1.VecBDTxzx = copy(WF2.VecBDTxzx)
  WF1.VecBDTxzz = copy(WF2.VecBDTxzz)
end
