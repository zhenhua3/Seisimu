function Onestep2D!(model::nspmod2d)
# For unsplit PML boundary simulation

  wf = model.wf
  fd = model.fd
  pml = model.pml

  vxbtxx = fd.Amvxbtxx.*(wf.txx*fd.Dxvxbtxx')
  vxbtxz = fd.Amvxbtxz.*(fd.Dzvxbtxz*wf.txz)
  vzbtzz = fd.Amvzbtzz.*(fd.Dzvzbtzz*wf.tzz)
  vzbtxz = fd.Amvzbtxz.*(wf.txz*fd.Dxvzbtxz')

  pml.PVxBTxx = (pml.bxVx .* pml.PVxBTxx + pml.PVxBTxx + pml.axVx .* vxbtxx)
  pml.PVxBTxz = (pml.bzVx .* pml.PVxBTxz + pml.PVxBTxz + pml.azVx .* vxbtxz)
  pml.PVzBTzz = (pml.bzVz .* pml.PVzBTzz + pml.PVzBTzz + pml.azVz .* vzbtzz)
  pml.PVzBTxz = (pml.bxVz .* pml.PVzBTxz + pml.PVzBTxz + pml.axVz .* vzbtxz)

  wf.vx = wf.vx + vxbtxx + vxbtxz + pml.PVxBTxx + pml.PVxBTxz
  wf.vz = wf.vz + vzbtzz + vzbtxz + pml.PVzBTzz + pml.PVzBTxz

  txxbvx = fd.Amtxxbvx.*(wf.vx*fd.Dxtxxbvx')
  txxbvz = fd.Amtxxbvz.*(fd.Dztxxbvz*wf.vz)
  tzzbvx = fd.Amtzzbvx.*(wf.vx*fd.Dxtzzbvx')
  tzzbvz = fd.Amtzzbvz.*(fd.Dztzzbvz*wf.vz)
  txzbvx = fd.Amtxzbvx.*(fd.Dztxzbvx*wf.vx)
  txzbvz = fd.Amtxzbvz.*(wf.vz*fd.Dxtxzbvz')

  pml.PTxxBVx = (pml.bxTxx .* pml.PTxxBVx + pml.PTxxBVx + pml.axTxx .* txxbvx)
  pml.PTxxBVz = (pml.bzTxx .* pml.PTxxBVz + pml.PTxxBVz + pml.azTxx .* txxbvz)
  pml.PTzzBVx = (pml.bxTzz .* pml.PTzzBVx + pml.PTzzBVx + pml.axTzz .* tzzbvx)
  pml.PTzzBVz = (pml.bzTzz .* pml.PTzzBVz + pml.PTzzBVz + pml.azTzz .* tzzbvz)
  pml.PTxzBVx = (pml.bzTxz .* pml.PTxzBVx + pml.PTxzBVx + pml.azTxz .* txzbvx)
  pml.PTxzBVz = (pml.bxTxz .* pml.PTxzBVz + pml.PTxzBVz + pml.axTxz .* txzbvz)

  wf.txx = wf.txx + txxbvx + txxbvz + pml.PTxxBVx + pml.PTxxBVz
  wf.tzz = wf.tzz + tzzbvx + tzzbvz + pml.PTzzBVx + pml.PTzzBVz
  wf.txz = wf.txz + txzbvx + txzbvz + pml.PTxzBVx + pml.PTxzBVz

end

function Onestep2D!(model::spmod2d)
  wf = model.wf
  fd = model.fd
  wf.vxx = fd.Amvxxbvxx.*wf.vxx + (wf.txxx+wf.txxz)*fd.Dxvxxbtxx'
  wf.vxz = fd.Amvxzbvxz.*wf.vxz + fd.Dzvxzbtxz*(wf.txzx+wf.txzz)
  wf.vzx = fd.Amvzxbvzx.*wf.vzx + (wf.txzx+wf.txzz)*fd.Dxvzxbtxz'
  wf.vzz = fd.Amvzzbvzz.*wf.vzz + fd.Dzvzzbtzz*(wf.tzzx+wf.tzzz)

  wf.txxx = fd.Amtxxxbtxxx.*wf.txxx + (wf.vxx+wf.vxz)*fd.Dxtxxxbvx'
  wf.txxz = fd.Amtxxzbtxxz.*wf.txxz + fd.Dztxxzbvz*(wf.vzx+wf.vzz)
  wf.tzzx = fd.Amtzzxbtzzx.*wf.tzzx + (wf.vxx+wf.vxz)*fd.Dxtzzxbvx'
  wf.tzzz = fd.Amtzzzbtzzz.*wf.tzzz + fd.Dztzzzbvz*(wf.vzx+wf.vzz)
  wf.txzx = fd.Amtxzxbtxzx.*wf.txzx + (wf.vzx+wf.vzz)*fd.Dxtxzxbvz'
  wf.txzz = fd.Amtxzzbtxzz.*wf.txzz + fd.Dztxzzbvx*(wf.vxx+wf.vxz)
end
