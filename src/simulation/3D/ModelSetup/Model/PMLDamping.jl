type dampCoef
    bxVx::SparseMatrixCSC{Float64, Int64}
    bzVx::SparseMatrixCSC{Float64, Int64}
    axVx::SparseMatrixCSC{Float64, Int64}
    azVx::SparseMatrixCSC{Float64, Int64}
    bxVz::SparseMatrixCSC{Float64, Int64}
    bzVz::SparseMatrixCSC{Float64, Int64}
    axVz::SparseMatrixCSC{Float64, Int64}
    azVz::SparseMatrixCSC{Float64, Int64}
    bxTxx::SparseMatrixCSC{Float64, Int64}
    bzTxx::SparseMatrixCSC{Float64, Int64}
    axTxx::SparseMatrixCSC{Float64, Int64}
    azTxx::SparseMatrixCSC{Float64, Int64}
    bxTzz::SparseMatrixCSC{Float64, Int64}
    bzTzz::SparseMatrixCSC{Float64, Int64}
    axTzz::SparseMatrixCSC{Float64, Int64}
    azTzz::SparseMatrixCSC{Float64, Int64}
    bxTxz::SparseMatrixCSC{Float64, Int64}
    bzTxz::SparseMatrixCSC{Float64, Int64}
    axTxz::SparseMatrixCSC{Float64, Int64}
    azTxz::SparseMatrixCSC{Float64, Int64}
    PVxBTxx::SparseMatrixCSC{Float64, Int64}
    PVxBTxz::SparseMatrixCSC{Float64, Int64}
    PVzBTzz::SparseMatrixCSC{Float64, Int64}
    PVzBTxz::SparseMatrixCSC{Float64, Int64}
    PTxxBVx::SparseMatrixCSC{Float64, Int64}
    PTxxBVz::SparseMatrixCSC{Float64, Int64}
    PTzzBVx::SparseMatrixCSC{Float64, Int64}
    PTzzBVz::SparseMatrixCSC{Float64, Int64}
    PTxzBVx::SparseMatrixCSC{Float64, Int64}
    PTxzBVz::SparseMatrixCSC{Float64, Int64}
end

function DampBound(WFSize::WF2DSize, medium::Model, ext::Int64, iflag::Int64)

    # This function provides PML damping boundaries
    # Outputs are damping matrixs of Vx, Vz, Txx, Tzz, Txz (2D)
    # Nz : Spatial grid number in Z direction
    # Nx : Spatial grid number in X direction
    # maxv : Maximum medium velocity
    # ext : PML boundary thinkness at one side (25 is used normally)
    # iflag : value is 1 for unlimited boundary or value is 2 for free surface boundary

    kx = 1
    kz = 1 # kx and kz dont have large effact to the results so we choose 1

    maxv = maximum(medium.VP)

    # Vx
    (nz,nx) = WFSize.N_Vx
    tmpdVxz,tmpdVxx,ax,az = DampCoeff(nz, nx, medium.NDep, medium.NHor, ext, maxv, iflag)
    (nzVxz,nxVxz) = size(tmpdVxz)
    (nzVxx,nxVxx) = size(tmpdVxx)
    bxVx = sparse(reshape(full(exp(-(tmpdVxx./kx.+ax).*medium.dt).-1),nzVxx*nxVxx))
    bzVx = sparse(reshape(full(exp(-(tmpdVxz./kz.+az).*medium.dt).-1),nzVxz*nxVxz))
    axVx = sparse(reshape(full(tmpdVxx./(kx.*(tmpdVxx.+kx.*ax))),nzVxx*nxVxx)).*(bxVx)
    azVx = sparse(reshape(full(tmpdVxz./(kz.*(tmpdVxz.+kz.*az))),nzVxz*nxVxz)).*(bzVx)
    PVxBTxx = spzeros(WFSize.N_BDVx[1]*WFSize.N_BDVx[2])
    PVxBTxz = spzeros(WFSize.N_BDVx[1]*WFSize.N_BDVx[2])

    # Vz
    (nz,nx) = WFSize.N_Vz
    tmpdVzz,tmpdVzx,ax,az = DampCoeff(nz, nx, medium.NDep, medium.NHor, ext, maxv, iflag)
    (nzVzz,nxVzz) = size(tmpdVzz)
    (nzVzx,nxVzx) = size(tmpdVzx)
    bxVz = sparse(reshape(full(exp(-(tmpdVzx/kx.+ax)*medium.dt).-1),nzVzx*nxVzx))
    bzVz = sparse(reshape(full(exp(-(tmpdVzz/kz.+az)*medium.dt).-1),nzVzz*nxVzz))
    axVz = sparse(reshape(full(tmpdVzx./(kx*(tmpdVzx+kx*ax))),nzVzx*nxVzx)).*(bxVz)
    azVz = sparse(reshape(full(tmpdVzz./(kz*(tmpdVzz+kz*az))),nzVzz*nxVzz)).*(bzVz)
    PVzBTzz = spzeros(WFSize.N_BDVz[1]*WFSize.N_BDVz[2])
    PVzBTxz = spzeros(WFSize.N_BDVz[1]*WFSize.N_BDVz[2])

    # Txx
    (nz,nx) = WFSize.N_Tpp
    tmpdTxxz,tmpdTxxx,ax,az = DampCoeff(nz, nx, medium.NDep, medium.NHor, ext, maxv, iflag)
    (nzTxxx,nxTxxx) = size(tmpdTxxx)
    (nzTxxz,nxTxxz) = size(tmpdTxxz)
    bxTxx = sparse(reshape(full(exp(-(tmpdTxxx/kx.+ax)*medium.dt).-1),nzTxxx*nxTxxx))
    bzTxx = sparse(reshape(full(exp(-(tmpdTxxz/kz.+az)*medium.dt).-1),nzTxxz*nxTxxz))
    axTxx = sparse(reshape(full(tmpdTxxx./(kx*(tmpdTxxx+kx*ax))),nzTxxx*nxTxxx)).*(bxTxx)
    azTxx = sparse(reshape(full(tmpdTxxz./(kz*(tmpdTxxz+kz*az))),nzTxxz*nxTxxz)).*(bzTxx)
    PTxxBVx = spzeros(WFSize.N_BDTpp[1]*WFSize.N_BDTpp[2])
    PTxxBVz = spzeros(WFSize.N_BDTpp[1]*WFSize.N_BDTpp[2])


    # Tzz
    (nz,nx) = WFSize.N_Tpp
    tmpdTzzz,tmpdTzzx,ax,az = DampCoeff(nz, nx, medium.NDep, medium.NHor, ext, maxv, iflag)
    (nzTzzx,nxTzzx) = size(tmpdTzzx)
    (nzTzzz,nxTzzz) = size(tmpdTzzz)
    bxTzz = sparse(reshape(full(exp(-(tmpdTzzx/kx.+ax)*medium.dt).-1),nzTzzx*nxTzzx))
    bzTzz = sparse(reshape(full(exp(-(tmpdTzzz/kz.+az)*medium.dt).-1),nzTzzz*nxTzzz))
    axTzz = sparse(reshape(full(tmpdTzzx./(kx*(tmpdTzzx+kx*ax))),nzTzzx*nxTzzx)).*(bxTzz)
    azTzz = sparse(reshape(full(tmpdTzzz./(kz*(tmpdTzzz+kz*az))),nzTzzz*nxTzzz)).*(bzTzz)
    PTzzBVx = spzeros(WFSize.N_BDTpp[1]*WFSize.N_BDTpp[2])
    PTzzBVz = spzeros(WFSize.N_BDTpp[1]*WFSize.N_BDTpp[2])

    # Txz
    (nz,nx) = WFSize.N_Txz
    tmpdTxzz,tmpdTxzx,ax,az = DampCoeff(nz, nx, medium.NDep, medium.NHor, ext, maxv, iflag)
    (nzTxzx,nxTxzx) = size(tmpdTxzx)
    (nzTxzz,nxTxzz) = size(tmpdTxzz)
    bxTxz = sparse(reshape(full(exp(-(tmpdTxzx/kx.+ax)*medium.dt).-1),nzTxzx*nxTxzx))
    bzTxz = sparse(reshape(full(exp(-(tmpdTxzz/kz.+az)*medium.dt).-1),nzTxzz*nxTxzz))
    axTxz = sparse(reshape(full(tmpdTxzx./(kx.*(tmpdTxzx.+kx*ax))),nzTxzx*nxTxzx)).*(bxTxz)
    azTxz = sparse(reshape(full(tmpdTxzz./(kz.*(tmpdTxzz.+kz*az))),nzTxzz*nxTxzz)).*(bzTxz)
    PTxzBVx = spzeros(WFSize.N_BDTxz[1]*WFSize.N_BDTxz[2])
    PTxzBVz = spzeros(WFSize.N_BDTxz[1]*WFSize.N_BDTxz[2])

    dpCoef = dampCoef(bxVx, bzVx, axVx, azVx,
                      bxVz, bzVz, axVz, azVz,
                      bxTxx, bzTxx, axTxx, azTxx,
                      bxTzz, bzTzz, axTzz, azTzz,
                      bxTxz, bzTxz, axTxz, azTxz,
                      PVxBTxx, PVxBTxz,
                      PVzBTzz, PVzBTxz,
                      PTxxBVx, PTxxBVz,
                      PTzzBVx, PTzzBVz,
                      PTxzBVx, PTxzBVz)

    return dpCoef
end

# This function introduce Damping Coefficient to each damping layer
function DampCoeff(nz::Int64, nx::Int64, Nz::Int64, Nx::Int64, ext::Int64, maxv::Float64, iflag::Int64)

    m = 0.25
    n = 0.75
    ll = 4

    # R is a parameter depend on the thickness of boundaries
    if ext == 5
       R = 0.01
    elseif ext == 10
       R = 0.001
    elseif ext == 20
       R = 0.0001
    else
       error("unsupported damping layers, choose from 5, 10 or 20")
    end

    a_max = 10*pi
    tmp_a = linspace(0,a_max,ext)

    if iflag == 2 #unlimited medium

        dWFx = spzeros(nz+2*ext,nx+2*ext)
        dWFz = spzeros(nz+2*ext,nx+2*ext)
        ax = spzeros(nz+2*ext,nx+2*ext)
        az = spzeros(nz+2*ext,nx+2*ext)

            if is(nx,Nx) == true

        for i = 1 : ext
            dWFx[:,i] = -maxv/ext * log(R) * (m*(ext+1-i )/ext + n*((ext+1-i )/ext)^ll)
            ax[:,i] = ax[:,i].+tmp_a[i]
        end

        for i = nx+ext+1 : nx+2*ext
            dWFx[:,i] = -maxv/ext * log(R) * (m*(i-nx-ext)/ext + n*((i-nx-ext)/ext)^ll)
            ax[:,i] = ax[:,i].+tmp_a[nx+2*ext-i+1]
        end

            elseif is(nx+1,Nx) == true

        for i = 1 : ext
            dWFx[:,i] = -maxv/ext * log(R) * (m*(ext+0.5-i)/ext + n*((ext+0.5-i)/ext)^ll)
            ax[:,i] = ax[:,i].+tmp_a[i]
        end
        for i = (nx+ext+1) : (nx+2*ext)
            dWFx[:,i] = -maxv/ext * log(R) * (m*(i-(nx)-ext-0.5)/ext + n*((i-(nx)-ext-0.5)/ext)^ll)
            ax[:,i] = ax[:,i].+tmp_a[nx+2*ext-i+1]
        end

      else error("please check your wavefield size")

            end

            if is(nz,Nz) == true

        for i = 1 : ext
            dWFz[i,:] = -maxv/ext * log(R) * (m*(ext+1-i )/ext + n*((ext+1-i )/ext)^ll)
            az[i,:] = az[i,:].+tmp_a[i]
        end
        for i = (nz+ext+1) : nz+2*ext
            dWFz[i,:] = -maxv/ext * log(R) * (m*(i-nz-ext)/ext + n*((i-nz-ext)/ext)^ll)
            az[i,:] = az[i,:].+tmp_a[nz+2*ext-i+1]
        end

            elseif is(nz+1,Nz) == true

        for i = 1 : ext
            dWFz[i,:] = -maxv/ext * log(R) * (m*(ext+0.5-i)/ext + n*((ext+0.5-i)/ext)^ll)
            az[i,:] = az[i,:].+tmp_a[i]
        end
        for i = nz+ext+1 : nz+2*ext
            dWFz[i,:] = -maxv/ext * log(R) * (m*(i-(nz)-ext-0.5)/ext + n*((i-(nz)-ext-0.5)/ext)^ll)
            az[i,:] = az[i,:].+tmp_a[nz+2*ext-i+1]
        end
        else error("please check your wavefield size")
            end

    elseif iflag == 1 # free surface

            dWFx = spzeros(nz+ext,nx+2*ext)
            dWFz = spzeros(nz+ext,nx+2*ext)
            ax = spzeros(nz+ext,nx+2*ext)
            az = spzeros(nz+ext,nx+2*ext)

            if is(nx,Nx) == true

        for i = 1 : ext
            dWFx[:,i] = -maxv/ext * log(R) * (m*(ext+1-i )/ext + n*((ext+1-i )/ext)^ll)
            ax[:,i] = ax[:,i].+tmp_a[i]
        end
        for i = nx+ext+1 : nx+2*ext
            dWFx[:,i] = -maxv/ext * log(R) * (m*(i-nx-ext)/ext + n*((i-nx-ext)/ext)^ll)
            ax[:,i] = ax[:,i].+tmp_a[nx+2*ext-i+1]
        end

            elseif is(nx+1,Nx) == true

        for i = 1 : ext
            dWFx[:,i] = -maxv/ext * log(R) * (m*(ext+0.5-i)/ext + n*((ext+0.5-i)/ext)^ll)
            ax[:,i] = ax[:,i].+tmp_a[i]
        end
        for i = (nx+ext+1) : (nx+2*ext)
            dWFx[:,i] = -maxv/ext * log(R) * (m*(i-(nx)-ext-0.5)/ext + n*((i-(nx)-ext-0.5)/ext)^ll)
            ax[:,i] = ax[:,i].+tmp_a[nx+2*ext-i+1]
        end
            else error("please check your wavefield size")
            end

            if is(nz,Nz) == true

        for i = (nz+1) : nz+ext
            dWFz[i,:] = -maxv/ext * log(R) * (m*(i-nz)/ext + n*((i-nz)/ext)^ll)
            az[i,:] = az[i,:].+tmp_a[nz+ext-i+1]
        end

            elseif is(nz+1,Nz) == true

        for i = nz+1 : nz+ext
            dWFz[i,:] = -maxv/ext * log(R) * (m*(i-(nz)-0.5)/ext + n*((i-(nz)-0.5)/ext)^ll)
            az[i,:] = az[i,:].+tmp_a[nz+ext-i+1]
        end
            else error("please check your wavefield size")
            end

        end

        return dWFz, dWFx, ax, az
end
