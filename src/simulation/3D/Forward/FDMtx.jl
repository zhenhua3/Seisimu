# ================= when apply split PML boudary condition ================ #
type SFDMtx
  VxxBTxx::SparseMatrixCSC{Float64, Int64}
  VxzBTxz::SparseMatrixCSC{Float64, Int64}
  VzzBTzz::SparseMatrixCSC{Float64, Int64}
  VzxBTxz::SparseMatrixCSC{Float64, Int64}
  TxxxBVx::SparseMatrixCSC{Float64, Int64}
  TxxzBVz::SparseMatrixCSC{Float64, Int64}
  TzzxBVx::SparseMatrixCSC{Float64, Int64}
  TzzzBVz::SparseMatrixCSC{Float64, Int64}
  TxzzBVx::SparseMatrixCSC{Float64, Int64}
  TxzxBVz::SparseMatrixCSC{Float64, Int64}
  TxxxBTxxx::SparseMatrixCSC{Float64, Int64}
  TxxzBTxxz::SparseMatrixCSC{Float64, Int64}
  TzzxBTzzx::SparseMatrixCSC{Float64, Int64}
  TzzzBTzzz::SparseMatrixCSC{Float64, Int64}
  TxzxBTxzx::SparseMatrixCSC{Float64, Int64}
  TxzzBTxzz::SparseMatrixCSC{Float64, Int64}
  VxxBVxx::SparseMatrixCSC{Float64, Int64}
  VxzBVxz::SparseMatrixCSC{Float64, Int64}
  VzxBVzx::SparseMatrixCSC{Float64, Int64}
  VzzBVzz::SparseMatrixCSC{Float64, Int64}
end

function Fdmtx(WFSize::WF2DSize, medium::Model, FDCoeff::AbstractArray, ext::Int64, iflag::Int64)

  maxv = maximum(medium.VP)

      # Vx #########################
      (nz,nx) = WFSize.N_Vx
      tmpdVxz,tmpdVxx,ax,az = DampCoeff(nz, nx, medium.NDep, medium.NHor, ext, maxv, iflag)
      tmpdVxz = tmpdVxz
      tmpdVxx = tmpdVxx
      AvgRhoVx = zeros(medium.PML_NHor,medium.PML_NHor-1)
        for i = 1:medium.PML_NHor-1
          AvgRhoVx[i:i+1,i] = AvgRhoVx[i:i+1,i] + [0.5,0.5]
        end
      rho_Vx = medium.Rho*AvgRhoVx

      # VxxBTxx
      (nz,nx) = WFSize.N_BDVx
      Am = zeros(nz*nx)
      for j = 1:nx
          for i = 1:nz
            Am[i+(j-1)*nz] = 1/(rho_Vx[i,j]/medium.dt + tmpdVxx[i,j]*rho_Vx[i,j]/2)/medium.dx
          end
      end
      (nz,nx) = WFSize.N_BDTpp
      DxTxx = Dx2D(nz, nx, medium.PML_NDep, medium.PML_NHor, FDCoeff)
      VxxBTxx = Am.*DxTxx

      #VxxBVxx
      (nz,nx) = WFSize.N_BDVx
      Am = zeros(nz*nx)
      for j = 1:nx
        for i = 1:nz
          a3 = rho_Vx[i,j]/medium.dt + tmpdVxx[i,j]*rho_Vx[i,j]/2
          Am[i+(j-1)*nz] = rho_Vx[i,j]/(medium.dt*a3) - tmpdVxx[i,j]*rho_Vx[i,j]/(2*a3)
        end
      end
      VxxBVxx = spdiagm(Am)

      #VxzBTxz
      (nz,nx) = WFSize.N_BDVx
      Am = zeros(nz*nx)
      for j = 1:nx
          for i = 1:nz
            Am[i+(j-1)*nz] = 1/(rho_Vx[i,j]/medium.dt + tmpdVxz[i,j]*rho_Vx[i,j]/2)/medium.dz
          end
      end
      (nz,nx) = WFSize.N_BDTxz
      DzTxz = Dz2D(nz, nx, medium.PML_NDep, medium.PML_NHor, FDCoeff)
      VxzBTxz = Am.*DzTxz

      #VxzBVxz
      (nz,nx) = WFSize.N_BDVx
      Am = zeros(nz*nx)
      for j = 1:nx
        for i = 1:nz
          a3 = rho_Vx[i,j]/medium.dt + tmpdVxz[i,j]*rho_Vx[i,j]/2
          Am[i+(j-1)*nz] = rho_Vx[i,j]/(medium.dt*a3) - tmpdVxz[i,j]*rho_Vx[i,j]/(2*a3)
        end
      end
      VxzBVxz = spdiagm(Am)

      #Vz ##########################
      (nz,nx) = WFSize.N_Vz
      tmpdVzz,tmpdVzx,ax,az = DampCoeff(nz, nx, medium.NDep, medium.NHor, ext, maxv, iflag)
      tmpdVzz = tmpdVzz
      tmpdVzx = tmpdVzx
      AvgRhoVz = zeros(medium.PML_NDep-1,medium.PML_NDep)
      for i = 1:medium.PML_NDep-1
        AvgRhoVz[i:i,i:i+1] = AvgRhoVz[i:i,i:i+1] .+ [0.5 0.5]
      end
      rho_Vz = AvgRhoVz*medium.Rho

      #VzzBTzz
      (nz,nx) = WFSize.N_BDVz
      Am = zeros(nz*nx)
      for j = 1:nx
          for i = 1:nz
            Am[i+(j-1)*nz] = 1/(rho_Vz[i,j]/medium.dt + tmpdVzz[i,j]*rho_Vz[i,j]/2)/medium.dz
          end
      end
      (nz,nx) = WFSize.N_BDTpp
      DzTzz = Dz2D(nz, nx, medium.PML_NDep, medium.PML_NHor, FDCoeff)
      VzzBTzz = Am.*DzTzz

      #VzzBVzz
      (nz,nx) = WFSize.N_BDVz
      Am = zeros(nz*nx)
      for j = 1:nx
        for i = 1:nz
          a3 = rho_Vz[i,j]/medium.dt + tmpdVzz[i,j]*rho_Vz[i,j]/2
          Am[i+(j-1)*nz] = rho_Vz[i,j]/(medium.dt*a3) - tmpdVzz[i,j]*rho_Vz[i,j]/(2*a3)
        end
      end
      VzzBVzz = spdiagm(Am)

      #VzxBTxz
      (nz,nx) = WFSize.N_BDVz
      Am = zeros(nz*nx)
      for j = 1:nx
          for i = 1:nz
            Am[i+(j-1)*nz] = 1/(rho_Vz[i,j]/medium.dt + tmpdVzx[i,j]*rho_Vz[i,j]/2)/medium.dx
          end
      end
      (nz,nx) = WFSize.N_BDTxz
      DxTxz = Dx2D(nz, nx, medium.PML_NDep, medium.PML_NHor, FDCoeff)
      VzxBTxz = Am.*DxTxz

      #VzxBVzx
      (nz,nx) = WFSize.N_BDVz
      Am = zeros(nz*nx)
      for j = 1:nx
        for i = 1:nz
          a3 = rho_Vz[i,j]/medium.dt + tmpdVzx[i,j]*rho_Vz[i,j]/2
          Am[i+(j-1)*nz] = rho_Vz[i,j]/(medium.dt*a3) - tmpdVzx[i,j]*rho_Vz[i,j]/(2*a3)
        end
      end
      VzxBVzx = spdiagm(Am)

      #Txx ##########################
      (nz,nx) = WFSize.N_Tpp
      tmpdTxxz,tmpdTxxx,ax,az = DampCoeff(nz, nx, medium.NDep, medium.NHor, ext, maxv, iflag)
      tmpdTxxz = tmpdTxxz
      tmpdTxxx = tmpdTxxx
      #TxxxBVx
      (nz,nx) = WFSize.N_BDTpp
      Am = zeros(nz*nx)
      for j = 1:nx
        for i = 1:nz
          Am[i+(j-1)*nz] = (medium.Lambda[i,j]+2*medium.Mu[i,j])/(1/medium.dt + tmpdTxxx[i,j]/2)/medium.dx
        end
      end
      (nz,nx) = WFSize.N_BDVx
      DxVx = Dx2D(nz, nx, medium.PML_NDep, medium.PML_NHor, FDCoeff)
      TxxxBVx = Am.*DxVx

      #TxxxBTxxx
      (nz,nx) = WFSize.N_BDTpp
      Am = zeros(nz*nx)
      for j = 1:nx
        for i = 1:nz
          a3 = 1/medium.dt + tmpdTxxx[i,j]/2
          Am[i+(j-1)*nz] = 1/(medium.dt*a3) - tmpdTxxx[i,j]/(2*a3)
        end
      end
      TxxxBTxxx = spdiagm(Am)

      #TxxzBVz
      (nz,nx) = WFSize.N_BDTpp
      Am = zeros(nz*nx)
      for j = 1:nx
        for i = 1:nz
          Am[i+(j-1)*nz] = (medium.Lambda[i,j])/(1/medium.dt + tmpdTxxx[i,j]/2)/medium.dz
        end
      end
      (nz,nx) = WFSize.N_BDVz
      DzVz = Dz2D(nz, nx, medium.PML_NDep, medium.PML_NHor, FDCoeff)
      TxxzBVz = Am.*DzVz

      #TxxzBTxxz
      (nz,nx) = WFSize.N_BDTpp
      Am = zeros(nz*nx)
      for j = 1:nx
        for i = 1:nz
          a3 = 1/medium.dt + tmpdTxxz[i,j]/2
          Am[i+(j-1)*nz] = 1/(medium.dt*a3) - tmpdTxxz[i,j]/(2*a3)
        end
      end
      TxxzBTxxz = spdiagm(Am)

      #Tzz ##########################
      (nz,nx) = WFSize.N_Tpp
      tmpdTzzz,tmpdTzzx,ax,az = DampCoeff(nz, nx, medium.NDep, medium.NHor, ext, maxv, iflag)
      tmpdTzzz = tmpdTzzz
      tmpdTzzx = tmpdTzzx
      #TzzzBVz
      (nz,nx) = WFSize.N_BDTpp
      Am = zeros(nz*nx)
      for j = 1:nx
        for i = 1:nz
          Am[i+(j-1)*nz] = (medium.Lambda[i,j]+2*medium.Mu[i,j])/(1/medium.dt + tmpdTzzz[i,j]/2)/medium.dz
        end
      end
      (nz,nx) = WFSize.N_BDVz
      DzVz = Dz2D(nz, nx, medium.PML_NDep, medium.PML_NHor, FDCoeff)
      TzzzBVz = Am.*DzVz

      #TzzzBTzzz
      (nz,nx) = WFSize.N_BDTpp
      Am = zeros(nz*nx)
      for j = 1:nx
        for i = 1:nz
          a3 = 1/medium.dt + tmpdTzzz[i,j]/2
          Am[i+(j-1)*nz] = 1/(medium.dt*a3) - tmpdTzzz[i,j]/(2*a3)
        end
      end
      TzzzBTzzz = spdiagm(Am)

      #TzzxBVx
      (nz,nx) = WFSize.N_BDTpp
      Am = zeros(nz*nx)
      for j = 1:nx
        for i = 1:nz
          Am[i+(j-1)*nz] = (medium.Lambda[i,j])/(1/medium.dt + tmpdTzzx[i,j]/2)/medium.dx
        end
      end
      (nz,nx) = WFSize.N_BDVx
      DxVx = Dx2D(nz, nx, medium.PML_NDep, medium.PML_NHor, FDCoeff)
      TzzxBVx = Am.*DxVx

      #TzzxBTzzx
      (nz,nx) = WFSize.N_BDTpp
      Am = zeros(nz*nx)
      for j = 1:nx
        for i = 1:nz
          a3 = 1/medium.dt + tmpdTzzx[i,j]/2
          Am[i+(j-1)*nz] = 1/(medium.dt*a3) - tmpdTzzx[i,j]/(2*a3)
        end
      end
      TzzxBTzzx = spdiagm(Am)

      #Txz ##########################
      (nz,nx) = WFSize.N_Txz
      tmpdTxzz,tmpdTxzx,ax,az = DampCoeff(nz, nx, medium.NDep, medium.NHor, ext, maxv, iflag)
      tmpdTxzz = tmpdTxzz
      tmpdTxzx = tmpdTxzx
      AvgMuTxzX = zeros(medium.PML_NHor,medium.PML_NHor-1)
      for i = 1:medium.PML_NHor-1
        AvgMuTxzX[i:i+1,i] = AvgMuTxzX[i:i+1,i] + [0.25,0.25]
      end
      AvgMuTxzZ = zeros(medium.PML_NDep-1,medium.PML_NDep)
      for i = 1:medium.PML_NDep-1
        AvgMuTxzZ[i:i,i:i+1] = AvgMuTxzZ[i:i,i:i+1] + [1 1]
      end
      mu_Txz = AvgMuTxzZ*medium.Mu*AvgMuTxzX

      #TxzzBVx
      (nz,nx) = WFSize.N_BDTxz
      Am = zeros(nz*nx)
      for j = 1:nx
        for i = 1:nz
          Am[i+(j-1)*nz] = (mu_Txz[i,j])/(1/medium.dt + tmpdTxzz[i,j]/2)/medium.dz
        end
      end
      (nz,nx) = WFSize.N_BDVx
      DzVx = Dz2D(nz, nx, medium.PML_NDep, medium.PML_NHor, FDCoeff)
      TxzzBVx = Am.*DzVx

      #TxzzBTxzz
      (nz,nx) = WFSize.N_BDTxz
      Am = zeros(nz*nx)
      for j = 1:nx
        for i = 1:nz
          a3 = 1/medium.dt + tmpdTxzz[i,j]/2
          Am[i+(j-1)*nz] = 1/(medium.dt*a3) - tmpdTxzz[i,j]/(2*a3)
        end
      end
      TxzzBTxzz = spdiagm(Am)

      #TxzxBVz
      (nz,nx) = WFSize.N_BDTxz
      Am = zeros(nz*nx)
      for j = 1:nx
        for i = 1:nz
          Am[i+(j-1)*nz] = (mu_Txz[i,j])/(1/medium.dt + tmpdTxzx[i,j]/2)/medium.dx
        end
      end
      (nz,nx) = WFSize.N_BDVz
      DxVz = Dx2D(nz, nx, medium.PML_NDep, medium.PML_NHor, FDCoeff)
      TxzxBVz = Am.*DxVz

      #TxzxBTxzx
      (nz,nx) = WFSize.N_BDTxz
      Am = zeros(nz*nx)
      for j = 1:nx
        for i = 1:nz
          a3 = 1/medium.dt + tmpdTxzx[i,j]/2
          Am[i+(j-1)*nz] = 1/(medium.dt*a3) - tmpdTxzx[i,j]/(2*a3)
        end
      end
      TxzxBTxzx = spdiagm(Am)

      fdmtx = SFDMtx(VxxBTxx, VxzBTxz, VzzBTzz, VzxBTxz, TxxxBVx, TxxzBVz, TzzxBVx, TzzzBVz, TxzzBVx, TxzxBVz, TxxxBTxxx, TxxzBTxxz, TzzxBTzzx, TzzzBTzzz, TxzxBTxzx, TxzzBTxzz, VxxBVxx, VxzBVxz, VzxBVzx, VzzBVzz)

end

# ================= when apply unsplit PML boudary condition ================ #
type NSFDMtx
VxBTxx::SparseMatrixCSC{Float64, Int64}
VxBTxz::SparseMatrixCSC{Float64, Int64}
VzBTzz::SparseMatrixCSC{Float64, Int64}
VzBTxz::SparseMatrixCSC{Float64, Int64}
TxxBVx::SparseMatrixCSC{Float64, Int64}
TxxBVz::SparseMatrixCSC{Float64, Int64}
TzzBVx::SparseMatrixCSC{Float64, Int64}
TzzBVz::SparseMatrixCSC{Float64, Int64}
TxzBVx::SparseMatrixCSC{Float64, Int64}
TxzBVz::SparseMatrixCSC{Float64, Int64}
end

function Fdmtx(WFSize::WF2DSize, medium::Model, FDCoeff::AbstractArray, ext::Int64)

    # if iflag == 2 #unlimited medium #################################################

    # Vx ##########################
    AvgRhoVx = zeros(medium.PML_NHor,medium.PML_NHor-1)
      for i = 1:medium.PML_NHor-1
        AvgRhoVx[i:i+1,i] = AvgRhoVx[i:i+1,i] + [0.5,0.5]
      end
    rho_Vx = medium.Rho*AvgRhoVx

    # VxBTxx
    (nz,nx) = WFSize.N_BDVx
    Am = zeros(nz*nx)
    for j = 1:nx
        for i = 1:nz
          Am[i+(j-1)*nz] = medium.dt./rho_Vx[i,j]/medium.dx
        end
    end
    (nz,nx) = WFSize.N_BDTpp
    DxTxx = Dx2D(nz, nx, medium.PML_NDep, medium.PML_NHor, FDCoeff)
    gc()
    VxBTxx = Am.*DxTxx

    #VxBTxz
    (nz,nx) = WFSize.N_BDVx
    Am = zeros(nz*nx)
    for j = 1:nx
        for i = 1:nz
          Am[i+(j-1)*nz] = medium.dt./rho_Vx[i,j]/medium.dz
        end
    end
    (nz,nx) = WFSize.N_BDTxz
    DzTxz = Dz2D(nz, nx, medium.PML_NDep, medium.PML_NHor, FDCoeff)
    gc()
    VxBTxz = Am.*DzTxz

    #Vz ##########################
    AvgRhoVz = zeros(medium.PML_NDep-1,medium.PML_NDep)
    for i = 1:medium.PML_NDep-1
      AvgRhoVz[i:i,i:i+1] = AvgRhoVz[i:i,i:i+1] .+ [0.5 0.5]
    end
    rho_Vz = AvgRhoVz*medium.Rho

    #VzBTzz
    (nz,nx) = WFSize.N_BDVz
    Am = zeros(nz*nx)
    for j = 1:nx
        for i = 1:nz
          Am[i+(j-1)*nz] = medium.dt./rho_Vz[i,j]/medium.dz
        end
    end
    (nz,nx) = WFSize.N_BDTpp
    DzTzz = Dz2D(nz, nx, medium.PML_NDep, medium.PML_NHor, FDCoeff)
    gc()
    VzBTzz = Am.*DzTzz

    #VzBTxz
    (nz,nx) = WFSize.N_BDVz
    Am = zeros(nz*nx)
    for j = 1:nx
        for i = 1:nz
          Am[i+(j-1)*nz] = medium.dt./rho_Vz[i,j]/medium.dx
        end
    end
    (nz,nx) = WFSize.N_BDTxz
    DxTxz = Dx2D(nz, nx, medium.PML_NDep, medium.PML_NHor, FDCoeff)
    gc()
    VzBTxz = Am.*DxTxz

    #Txx ##########################

    #TxxBVx
    (nz,nx) = WFSize.N_BDTpp
    Am = zeros(nz*nx)
    for j = 1:nx
      for i = 1:nz
        Am[i+(j-1)*nz] = medium.dt.*(medium.Lambda[i,j]+2*medium.Mu[i,j])/medium.dx
      end
    end
    (nz,nx) = WFSize.N_BDVx
    DxVx = Dx2D(nz, nx, medium.PML_NDep, medium.PML_NHor, FDCoeff)
    gc()
    TxxBVx = Am.*DxVx

    #TxxBVz
    (nz,nx) = WFSize.N_BDTpp
    Am = zeros(nz*nx)
    for j = 1:nx
      for i = 1:nz
        Am[i+(j-1)*nz] = medium.dt.*(medium.Lambda[i,j])/medium.dz
      end
    end
    (nz,nx) = WFSize.N_BDVz
    DzVz = Dz2D(nz, nx, medium.PML_NDep, medium.PML_NHor, FDCoeff)
    gc()
    TxxBVz = Am.*DzVz

    #Tzz ##########################

    #TzzBVz
    (nz,nx) = WFSize.N_BDTpp
    Am = zeros(nz*nx)
    for j = 1:nx
      for i = 1:nz
        Am[i+(j-1)*nz] = medium.dt.*(medium.Lambda[i,j]+2*medium.Mu[i,j])/medium.dz
      end
    end
    (nz,nx) = WFSize.N_BDVz
    DzVz = Dz2D(nz, nx, medium.PML_NDep, medium.PML_NHor, FDCoeff)
    gc()
    TzzBVz = Am.*DzVz

    #TzzBVx
    (nz,nx) = WFSize.N_BDTpp
    Am = zeros(nz*nx)
    for j = 1:nx
      for i = 1:nz
        Am[i+(j-1)*nz] = medium.dt.*(medium.Lambda[i,j])/medium.dx
      end
    end
    (nz,nx) = WFSize.N_BDVx
    DxVx = Dx2D(nz, nx, medium.PML_NDep, medium.PML_NHor, FDCoeff)
    gc()
    TzzBVx = Am.*DxVx

    #Txz ##########################
    AvgMuTxzX = zeros(medium.PML_NHor,medium.PML_NHor-1)
    for i = 1:medium.PML_NHor-1
      AvgMuTxzX[i:i+1,i] = AvgMuTxzX[i:i+1,i] + [0.25,0.25]
    end
    AvgMuTxzZ = zeros(medium.PML_NDep-1,medium.PML_NDep)
    for i = 1:medium.PML_NDep-1
      AvgMuTxzZ[i:i,i:i+1] = AvgMuTxzZ[i:i,i:i+1] + [1 1]
    end
    mu_Txz = AvgMuTxzZ*medium.Mu*AvgMuTxzX

    #TxzBVx
    (nz,nx) = WFSize.N_BDTxz
    Am = zeros(nz*nx)
    for j = 1:nx
      for i = 1:nz
        Am[i+(j-1)*nz] = medium.dt.*(mu_Txz[i,j])/medium.dz
      end
    end
    (nz,nx) = WFSize.N_BDVx
    DzVx = Dz2D(nz, nx, medium.PML_NDep, medium.PML_NHor, FDCoeff)
    gc()
    TxzBVx = Am.*DzVx

    #TxzBVz
    (nz,nx) = WFSize.N_BDTxz
    Am = zeros(nz*nx)
    for j = 1:nx
      for i = 1:nz
        Am[i+(j-1)*nz] = medium.dt.*(mu_Txz[i,j])/medium.dx
      end
    end
    (nz,nx) = WFSize.N_BDVz
    DxVz = Dx2D(nz, nx, medium.PML_NDep, medium.PML_NHor, FDCoeff)
    gc()
    TxzBVz = Am.*DxVz

fdmtx = NSFDMtx(VxBTxx, VxBTxz, VzBTzz, VzBTxz, TxxBVx, TxxBVz, TzzBVx, TzzBVz, TxzBVx, TxzBVz)
end
