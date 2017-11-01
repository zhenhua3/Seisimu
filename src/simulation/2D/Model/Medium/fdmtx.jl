# ================= when apply split PML boudary condition ================ #
type spfdmtx2d
   Dxvxxbtxx::SparseMatrixCSC{Float64, Int64}
   Amvxxbtxx::Array{Float64,2}
   Dzvxzbtxz::SparseMatrixCSC{Float64, Int64}
   Amvxzbtxz::Array{Float64,2}
   Dzvzzbtzz::SparseMatrixCSC{Float64, Int64}
   Amvzzbtzz::Array{Float64,2}
   Dxvzxbtxz::SparseMatrixCSC{Float64, Int64}
   Amvzxbtxz::Array{Float64,2}
   Dxtxxxbvx::SparseMatrixCSC{Float64, Int64}
   Amtxxxbvx::Array{Float64,2}
   Dztxxzbvz::SparseMatrixCSC{Float64, Int64}
   Amtxxzbvz::Array{Float64,2}
   Dxtzzxbvx::SparseMatrixCSC{Float64, Int64}
   Amtzzxbvx::Array{Float64,2}
   Dztzzzbvz::SparseMatrixCSC{Float64, Int64}
   Amtzzzbvz::Array{Float64,2}
   Dztxzzbvx::SparseMatrixCSC{Float64, Int64}
   Amtxzzbvx::Array{Float64,2}
   Dxtxzxbvz::SparseMatrixCSC{Float64, Int64}
   Amtxzxbvz::Array{Float64,2}
   Amtxxxbtxxx::SparseMatrixCSC{Float64, Int64}
   Amtxxzbtxxz::SparseMatrixCSC{Float64, Int64}
   Amtzzxbtzzx::SparseMatrixCSC{Float64, Int64}
   Amtzzzbtzzz::SparseMatrixCSC{Float64, Int64}
   Amtxzxbtxzx::SparseMatrixCSC{Float64, Int64}
   Amtxzzbtxzz::SparseMatrixCSC{Float64, Int64}
   Amvxxbvxx::SparseMatrixCSC{Float64, Int64}
   Amvxzbvxz::SparseMatrixCSC{Float64, Int64}
   Amvzxbvzx::SparseMatrixCSC{Float64, Int64}
   Amvzzbvzz::SparseMatrixCSC{Float64, Int64}
end

function fdmtx(nwf::nwf2d, medium::medium2d, FDCoeff::AbstractArray, ext::Int64, iflag::Int64)

  maxv = maximum(medium.pvel)

      # Vx #########################
      (nz,nx) = nwf.nvx
      tmpdVxz,tmpdVxx,ax,az = DampCoeff(nz, nx, medium.nDZ, medium.nHX, ext, maxv, iflag)
      tmpdVxz = tmpdVxz
      tmpdVxx = tmpdVxx
      AvgRhoVx = zeros(medium.BDnHX,medium.BDnHX-1)
        for i = 1:medium.BDnHX-1
          AvgRhoVx[i:i+1,i] = AvgRhoVx[i:i+1,i] + [0.5,0.5]
        end
      rho_Vx = medium.rho*AvgRhoVx

      # VxxBTxx
      (nz,nx) = nwf.BDnvx
      AmVxxBTxx = zeros(nz,nx)
      for j = 1:nx
          for i = 1:nz
            AmVxxBTxx[i,j] = 1/(rho_Vx[i,j]/medium.dt + tmpdVxx[i,j]*rho_Vx[i,j]/2)/medium.dx
          end
      end
      (nz,nx) = nwf.BDntpp
      DxVxxBTxx = Dx2D(nz, nx, medium.BDnDZ, medium.BDnHX, FDCoeff)

      #VxxBVxx
      (nz,nx) = nwf.BDnvx
      AmVxxBVxx = zeros(nz,nx)
      for j = 1:nx
        for i = 1:nz
          a3 = rho_Vx[i,j]/medium.dt + tmpdVxx[i,j]*rho_Vx[i,j]/2
          AmVxxBVxx[i,j] = rho_Vx[i,j]/(medium.dt*a3) - tmpdVxx[i,j]*rho_Vx[i,j]/(2*a3)
        end
      end

      #VxzBTxz
      (nz,nx) = nwf.BDnvx
      AmVxzBTxz = zeros(nz,nx)
      for j = 1:nx
          for i = 1:nz
            AmVxzBTxz[i,j] = 1/(rho_Vx[i,j]/medium.dt + tmpdVxz[i,j]*rho_Vx[i,j]/2)/medium.dz
          end
      end
      (nz,nx) = nwf.BDntxz
      DzVxzBTxz = Dz2D(nz, nx, medium.BDnDZ, medium.BDnHX, FDCoeff)

      #VxzBVxz
      (nz,nx) = nwf.BDnvx
      AmVxzBVxz = zeros(nz,nx)
      for j = 1:nx
        for i = 1:nz
          a3 = rho_Vx[i,j]/medium.dt + tmpdVxz[i,j]*rho_Vx[i,j]/2
          AmVxzBVxz[i,j] = rho_Vx[i,j]/(medium.dt*a3) - tmpdVxz[i,j]*rho_Vx[i,j]/(2*a3)
        end
      end

      #Vz ##########################
      (nz,nx) = nwf.nvz
      tmpdVzz,tmpdVzx,ax,az = DampCoeff(nz, nx, medium.nDZ, medium.nHX, ext, maxv, iflag)
      tmpdVzz = tmpdVzz
      tmpdVzx = tmpdVzx
      AvgRhoVz = zeros(medium.BDnDZ-1,medium.BDnDZ)
      for i = 1:medium.BDnDZ-1
        AvgRhoVz[i:i,i:i+1] = AvgRhoVz[i:i,i:i+1] .+ [0.5 0.5]
      end
      rho_Vz = AvgRhoVz*medium.rho

      #VzzBTzz
      (nz,nx) = nwf.BDnvz
      AmVzzBTzz = zeros(nz,nx)
      for j = 1:nx
          for i = 1:nz
            AmVzzBTzz[i,j] = 1/(rho_Vz[i,j]/medium.dt + tmpdVzz[i,j]*rho_Vz[i,j]/2)/medium.dz
          end
      end
      (nz,nx) = nwf.BDntpp
      DzVzzBTzz = Dz2D(nz, nx, medium.BDnDZ, medium.BDnHX, FDCoeff)

      #VzzBVzz
      (nz,nx) = nwf.BDnvz
      AmVzzBVzz = zeros(nz,nx)
      for j = 1:nx
        for i = 1:nz
          a3 = rho_Vz[i,j]/medium.dt + tmpdVzz[i,j]*rho_Vz[i,j]/2
          AmVzzBVzz[i+(j-1)*nz] = rho_Vz[i,j]/(medium.dt*a3) - tmpdVzz[i,j]*rho_Vz[i,j]/(2*a3)
        end
      end

      #VzxBTxz
      (nz,nx) = nwf.BDnvz
      AmVzxBTxz = zeros(nz,nx)
      for j = 1:nx
          for i = 1:nz
            AmVzxBTxz[i,j] = 1/(rho_Vz[i,j]/medium.dt + tmpdVzx[i,j]*rho_Vz[i,j]/2)/medium.dx
          end
      end
      (nz,nx) = nwf.BDntxz
      DxVzxBTxz = Dx2D(nz, nx, medium.BDnDZ, medium.BDnHX, FDCoeff)

      #VzxBVzx
      (nz,nx) = nwf.BDnvz
      AmVzxBVzx = zeros(nz,nx)
      for j = 1:nx
        for i = 1:nz
          a3 = rho_Vz[i,j]/medium.dt + tmpdVzx[i,j]*rho_Vz[i,j]/2
          AmVzxBVzx[i,j] = rho_Vz[i,j]/(medium.dt*a3) - tmpdVzx[i,j]*rho_Vz[i,j]/(2*a3)
        end
      end

      #Txx ##########################
      (nz,nx) = nwf.ntpp
      tmpdTxxz,tmpdTxxx,ax,az = DampCoeff(nz, nx, medium.nDZ, medium.nHX, ext, maxv, iflag)
      tmpdTxxz = tmpdTxxz
      tmpdTxxx = tmpdTxxx
      #TxxxBVx
      (nz,nx) = nwf.BDntpp
      AmTxxxBVx = zeros(nz,nx)
      for j = 1:nx
        for i = 1:nz
          AmTxxxBVx[i,j] = (medium.lambda[i,j]+2*medium.mu[i,j])/(1/medium.dt + tmpdTxxx[i,j]/2)/medium.dx
        end
      end
      (nz,nx) = nwf.BDnvx
      DxTxxxBVx = Dx2D(nz, nx, medium.BDnDZ, medium.BDnHX, FDCoeff)

      #TxxxBTxxx
      (nz,nx) = nwf.BDntpp
      AmTxxxBTxxx = zeros(nz,nx)
      for j = 1:nx
        for i = 1:nz
          a3 = 1/medium.dt + tmpdTxxx[i,j]/2
          AmTxxxBTxxx[i,j] = 1/(medium.dt*a3) - tmpdTxxx[i,j]/(2*a3)
        end
      end

      #TxxzBVz
      (nz,nx) = nwf.BDntpp
      AmTxxzBVz = zeros(nz,nx)
      for j = 1:nx
        for i = 1:nz
          AmTxxzBVz[i,j] = (medium.lambda[i,j])/(1/medium.dt + tmpdTxxx[i,j]/2)/medium.dz
        end
      end
      (nz,nx) = nwf.BDnvz
      DzTxxzBVz = Dz2D(nz, nx, medium.BDnDZ, medium.BDnHX, FDCoeff)

      #TxxzBTxxz
      (nz,nx) = nwf.BDntpp
      AmTxxzBTxxz = zeros(nz,nx)
      for j = 1:nx
        for i = 1:nz
          a3 = 1/medium.dt + tmpdTxxz[i,j]/2
          AmTxxzBTxxz[i,j] = 1/(medium.dt*a3) - tmpdTxxz[i,j]/(2*a3)
        end
      end

      #Tzz ##########################
      (nz,nx) = nwf.ntpp
      tmpdTzzz,tmpdTzzx,ax,az = DampCoeff(nz, nx, medium.nDZ, medium.nHX, ext, maxv, iflag)
      tmpdTzzz = tmpdTzzz
      tmpdTzzx = tmpdTzzx
      #TzzzBVz
      (nz,nx) = nwf.BDntpp
      AmTzzzBVz = zeros(nz,nx)
      for j = 1:nx
        for i = 1:nz
          AmTzzzBVz[i,j] = (medium.lambda[i,j]+2*medium.mu[i,j])/(1/medium.dt + tmpdTzzz[i,j]/2)/medium.dz
        end
      end
      (nz,nx) = nwf.BDnvz
      DzTzzzBVz = Dz2D(nz, nx, medium.BDnDZ, medium.BDnHX, FDCoeff)

      #TzzzBTzzz
      (nz,nx) = nwf.BDntpp
      AmTzzzBTzzz = zeros(nz,nx)
      for j = 1:nx
        for i = 1:nz
          a3 = 1/medium.dt + tmpdTzzz[i,j]/2
          AmTzzzBTzzz[i,j] = 1/(medium.dt*a3) - tmpdTzzz[i,j]/(2*a3)
        end
      end

      #TzzxBVx
      (nz,nx) = nwf.BDntpp
      AmTzzxBVx = zeros(nz,nx)
      for j = 1:nx
        for i = 1:nz
          AmTzzxBVx[i,j] = (medium.lambda[i,j])/(1/medium.dt + tmpdTzzx[i,j]/2)/medium.dx
        end
      end
      (nz,nx) = nwf.BDnvx
      DxTzzxBVx = Dx2D(nz, nx, medium.BDnDZ, medium.BDnHX, FDCoeff)

      #TzzxBTzzx
      (nz,nx) = nwf.BDntpp
      AmTzzxBTzzx = zeros(nz,nx)
      for j = 1:nx
        for i = 1:nz
          a3 = 1/medium.dt + tmpdTzzx[i,j]/2
          AmTzzxBTzzx[i,j] = 1/(medium.dt*a3) - tmpdTzzx[i,j]/(2*a3)
        end
      end

      #Txz ##########################
      (nz,nx) = nwf.ntxz
      tmpdTxzz,tmpdTxzx,ax,az = DampCoeff(nz, nx, medium.nDZ, medium.nHX, ext, maxv, iflag)
      tmpdTxzz = tmpdTxzz
      tmpdTxzx = tmpdTxzx
      AvgMuTxzX = zeros(medium.BDnHX,medium.BDnHX-1)
      for i = 1:medium.BDnHX-1
        AvgMuTxzX[i:i+1,i] = AvgMuTxzX[i:i+1,i] + [0.25,0.25]
      end
      AvgMuTxzZ = zeros(medium.BDnDZ-1,medium.BDnDZ)
      for i = 1:medium.BDnDZ-1
        AvgMuTxzZ[i:i,i:i+1] = AvgMuTxzZ[i:i,i:i+1] + [1 1]
      end
      mu_Txz = AvgMuTxzZ*medium.mu*AvgMuTxzX

      #TxzzBVx
      (nz,nx) = nwf.BDntxz
      Am = zeros(nz,nx)
      for j = 1:nx
        for i = 1:nz
          AmTxzzBVx[i,j] = (mu_Txz[i,j])/(1/medium.dt + tmpdTxzz[i,j]/2)/medium.dz
        end
      end
      (nz,nx) = nwf.BDnvx
      DzTxzzBVx = Dz2D(nz, nx, medium.BDnDZ, medium.BDnHX, FDCoeff)

      #TxzzBTxzz
      (nz,nx) = nwf.BDntxz
      AmTxzzBTxzz = zeros(nz,nx)
      for j = 1:nx
        for i = 1:nz
          a3 = 1/medium.dt + tmpdTxzz[i,j]/2
          AmTxzzBTxzz[i,j] = 1/(medium.dt*a3) - tmpdTxzz[i,j]/(2*a3)
        end
      end

      #TxzxBVz
      (nz,nx) = nwf.BDntxz
      AmTxzxBVz = zeros(nz,nx)
      for j = 1:nx
        for i = 1:nz
          AmTxzxBVz[i,j] = (mu_Txz[i,j])/(1/medium.dt + tmpdTxzx[i,j]/2)/medium.dx
        end
      end
      (nz,nx) = nwf.BDnvz
      DxTxzxBVz = Dx2D(nz, nx, medium.BDnDZ, medium.BDnHX, FDCoeff)

      #TxzxBTxzx
      (nz,nx) = nwf.BDntxz
      AmTxzxBTxzx = zeros(nz,nx)
      for j = 1:nx
        for i = 1:nz
          a3 = 1/medium.dt + tmpdTxzx[i,j]/2
          AmTxzxBTxzx[i,j] = 1/(medium.dt*a3) - tmpdTxzx[i,j]/(2*a3)
        end
      end

      return spfdmtx2d(
      DxVxxBTxx, AmVxxBTxx, DzVxzBTxz, AmVxzBTxz, DzVzzBTzz,
      AmVzzBTzz, DxVzxBTxz, AmVzxBTxz,
      DxTxxxBVx, AmTxxxBVx, DzTxxzBVz, AmTxxzBVz,
      DxTzzxBVx, AmTzzxBVx, DzTzzzBVz, AmTzzzBVz,
      DzTxzzBVx, AmTxzzBVx, DxTxzxBVz, AmTxzxBVz,
      AmTxxxBTxxx, AmTxxzBTxxz,
      AmTzzxBTzzx, AmTzzzBTzzz,
      AmTxzxBTxzx, AmTxzzBTxzz,
      AmVxxBVxx, AmVxzBVxz,
      AmVzxBVzx, AmVzzBVzz)

end






# ================= when apply unsplit PML boudary condition ================ #
type nspfdmtx2d
    Dxvxbtxx::SparseMatrixCSC{Float64, Int64}
    Amvxbtxx::Array{Float64,2}
    Dzvxbtxz::SparseMatrixCSC{Float64, Int64}
    Amvxbtxz::Array{Float64,2}
    Dzvzbtzz::SparseMatrixCSC{Float64, Int64}
    Amvzbtzz::Array{Float64,2}
    Dxvzbtxz::SparseMatrixCSC{Float64, Int64}
    Amvzbtxz::Array{Float64,2}
    Dxtxxbvx::SparseMatrixCSC{Float64, Int64}
    Amtxxbvx::Array{Float64,2}
    Dztxxbvz::SparseMatrixCSC{Float64, Int64}
    Amtxxbvz::Array{Float64,2}
    Dxtzzbvx::SparseMatrixCSC{Float64, Int64}
    Amtzzbvx::Array{Float64,2}
    Dztzzbvz::SparseMatrixCSC{Float64, Int64}
    Amtzzbvz::Array{Float64,2}
    Dztxzbvx::SparseMatrixCSC{Float64, Int64}
    Amtxzbvx::Array{Float64,2}
    Dxtxzbvz::SparseMatrixCSC{Float64, Int64}
    Amtxzbvz::Array{Float64,2}
end

function fdmtx(nwf::nwf2d, medium::medium2d, FDCoeff::AbstractArray, ext::Int64)

    # if iflag == 2 #unlimited medium #################################################

    # Vx ##########################
    AvgRhoVx = zeros(medium.BDnHX,medium.BDnHX-1)
      for i = 1:medium.BDnHX-1
        AvgRhoVx[i:i+1,i] = AvgRhoVx[i:i+1,i] + [0.5,0.5]
      end
    rho_Vx = medium.rho*AvgRhoVx

    # VxBTxx
    (nz,nx) = nwf.BDnvx
    AmVxBTxx = zeros(nz,nx)
    for j = 1:nx
        for i = 1:nz
          AmVxBTxx[i,j] = medium.dt./rho_Vx[i,j]/medium.dx
        end
    end
    (nz,nx) = nwf.BDntpp
    DxVxBTxx = Dx2D(nz, nx, medium.BDnDZ, medium.BDnHX, FDCoeff)

    #VxBTxz
    (nz,nx) = nwf.BDnvx
    AmVxBTxz = zeros(nz,nx)
    for j = 1:nx
        for i = 1:nz
          AmVxBTxz[i,j] = medium.dt./rho_Vx[i,j]/medium.dz
        end
    end
    (nz,nx) = nwf.BDntxz
    DzVxBTxz = Dz2D(nz, nx, medium.BDnDZ, medium.BDnHX, FDCoeff)

    #Vz ##########################
    AvgRhoVz = zeros(medium.BDnDZ-1,medium.BDnDZ)
    for i = 1:medium.BDnDZ-1
      AvgRhoVz[i:i,i:i+1] = AvgRhoVz[i:i,i:i+1] .+ [0.5 0.5]
    end
    rho_Vz = AvgRhoVz*medium.rho

    #VzBTzz
    (nz,nx) = nwf.BDnvz
    AmVzBTzz = zeros(nz,nx)
    for j = 1:nx
        for i = 1:nz
          AmVzBTzz[i,j] = medium.dt./rho_Vz[i,j]/medium.dz
        end
    end
    (nz,nx) = nwf.BDntpp
    DzVzBTzz = Dz2D(nz, nx, medium.BDnDZ, medium.BDnHX, FDCoeff)

    #VzBTxz
    (nz,nx) = nwf.BDnvz
    AmVzBTxz = zeros(nz,nx)
    for j = 1:nx
        for i = 1:nz
          AmVzBTxz[i,j] = medium.dt./rho_Vz[i,j]/medium.dx
        end
    end
    (nz,nx) = nwf.BDntxz
    DxVzBTxz = Dx2D(nz, nx, medium.BDnDZ, medium.BDnHX, FDCoeff)

    #Txx ##########################

    #TxxBVx
    (nz,nx) = nwf.BDntpp
    AmTxxBVx = zeros(nz,nx)
    for j = 1:nx
      for i = 1:nz
        AmTxxBVx[i,j] = medium.dt.*(medium.lambda[i,j]+2*medium.mu[i,j])/medium.dx
      end
    end
    (nz,nx) = nwf.BDnvx
    DxTxxBVx = Dx2D(nz, nx, medium.BDnDZ, medium.BDnHX, FDCoeff)

    #TxxBVz
    (nz,nx) = nwf.BDntpp
    AmTxxBVz = zeros(nz,nx)
    for j = 1:nx
      for i = 1:nz
        AmTxxBVz[i,j] = medium.dt.*(medium.lambda[i,j])/medium.dz
      end
    end
    (nz,nx) = nwf.BDnvz
    DzTxxBVz = Dz2D(nz, nx, medium.BDnDZ, medium.BDnHX, FDCoeff)

    #Tzz ##########################

    #TzzBVz
    (nz,nx) = nwf.BDntpp
    AmTzzBVz = zeros(nz,nx)
    for j = 1:nx
      for i = 1:nz
        AmTzzBVz[i,j] = medium.dt.*(medium.lambda[i,j]+2*medium.mu[i,j])/medium.dz
      end
    end
    (nz,nx) = nwf.BDnvz
    DzTzzBVz = Dz2D(nz, nx, medium.BDnDZ, medium.BDnHX, FDCoeff)

    #TzzBVx
    (nz,nx) = nwf.BDntpp
    AmTzzBVx = zeros(nz,nx)
    for j = 1:nx
      for i = 1:nz
        AmTzzBVx[i,j] = medium.dt.*(medium.lambda[i,j])/medium.dx
      end
    end
    (nz,nx) = nwf.BDnvx
    DxTzzBVx = Dx2D(nz, nx, medium.BDnDZ, medium.BDnHX, FDCoeff)

    #Txz ##########################
    AvgMuTxzX = zeros(medium.BDnHX,medium.BDnHX-1)
    for i = 1:medium.BDnHX-1
      AvgMuTxzX[i:i+1,i] = AvgMuTxzX[i:i+1,i] + [0.25,0.25]
    end
    AvgMuTxzZ = zeros(medium.BDnDZ-1,medium.BDnDZ)
    for i = 1:medium.BDnDZ-1
      AvgMuTxzZ[i:i,i:i+1] = AvgMuTxzZ[i:i,i:i+1] + [1 1]
    end
    mu_Txz = AvgMuTxzZ*medium.mu*AvgMuTxzX

    #TxzBVx
    (nz,nx) = nwf.BDntxz
    AmTxzBVx = zeros(nz,nx)
    for j = 1:nx
      for i = 1:nz
        AmTxzBVx[i,j] = medium.dt.*(mu_Txz[i,j])/medium.dz
      end
    end
    (nz,nx) = nwf.BDnvx
    DzTxzBVx = Dz2D(nz, nx, medium.BDnDZ, medium.BDnHX, FDCoeff)

    #TxzBVz
    (nz,nx) = nwf.BDntxz
    AmTxzBVz = zeros(nz,nx)
    for j = 1:nx
      for i = 1:nz
        AmTxzBVz[i,j] = medium.dt.*(mu_Txz[i,j])/medium.dx
      end
    end
    (nz,nx) = nwf.BDnvz
    DxTxzBVz = Dx2D(nz, nx, medium.BDnDZ, medium.BDnHX, FDCoeff)

    return nspfdmtx2d(
    DxVxBTxx, AmVxBTxx, DzVxBTxz, AmVxBTxz, DzVzBTzz, AmVzBTzz, DxVzBTxz, AmVzBTxz,
    DxTxxBVx, AmTxxBVx, DzTxxBVz, AmTxxBVz, DxTzzBVx, AmTzzBVx, DzTzzBVz, AmTzzBVz,
    DzTxzBVx, AmTxzBVx, DxTxzBVz, AmTxzBVz)
end
