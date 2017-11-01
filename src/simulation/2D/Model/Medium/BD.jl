type BD2d
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

function BD(nwf::nwf2d, medium::medium2d, ext::Int64, iflag::Int64)

    # This function provides PML damping boundaries
    # Outputs are damping matrixs of Vx, Vz, Txx, Tzz, Txz (2D)
    # nDZ : Spatial grid number in Z direction
    # nHX : Spatial grid number in X direction
    # maxv : Maximum medium velocity
    # ext : PML boundary thinkness at one side (25 is used normally)
    # iflag : value is 1 for unlimited boundary or value is 2 for free surface boundary

    kx = 1
    kz = 1 # kx and kz dont have large effact to the results so we choose 1

    maxv = maximum(medium.pvel)

    # Vx
    (nz,nx) = nwf.nvx
    tmpdVxz,tmpdVxx,ax,az = DampCoeff(nz, nx, medium.nDZ, medium.nHX, ext, maxv, iflag)
    (nzVxz,nxVxz) = size(tmpdVxz)
    (nzVxx,nxVxx) = size(tmpdVxx)
    bxVx = exp.(-(tmpdVxx./kx.+ax).*medium.dt).-1
    bzVx = exp.(-(tmpdVxz./kz.+az).*medium.dt).-1
    axVx = tmpdVxx./(kx.*(tmpdVxx.+kx.*ax)).*bxVx
    azVx = tmpdVxz./(kz.*(tmpdVxz.+kz.*az)).*bzVx
    PVxBTxx = spzeros(nwf.BDnvx[1],nwf.BDnvx[2])
    PVxBTxz = spzeros(nwf.BDnvx[1],nwf.BDnvx[2])

    # Vz
    (nz,nx) = nwf.nvz
    tmpdVzz,tmpdVzx,ax,az = DampCoeff(nz, nx, medium.nDZ, medium.nHX, ext, maxv, iflag)
    (nzVzz,nxVzz) = size(tmpdVzz)
    (nzVzx,nxVzx) = size(tmpdVzx)
    bxVz = exp.(-(tmpdVzx/kx.+ax)*medium.dt).-1
    bzVz = exp.(-(tmpdVzz/kz.+az)*medium.dt).-1
    axVz = tmpdVzx./(kx*(tmpdVzx+kx*ax)).*bxVz
    azVz = tmpdVzz./(kz*(tmpdVzz+kz*az)).*bzVz
    PVzBTzz = spzeros(nwf.BDnvz[1],nwf.BDnvz[2])
    PVzBTxz = spzeros(nwf.BDnvz[1],nwf.BDnvz[2])

    # Txx
    (nz,nx) = nwf.ntpp
    tmpdTxxz,tmpdTxxx,ax,az = DampCoeff(nz, nx, medium.nDZ, medium.nHX, ext, maxv, iflag)
    (nzTxxx,nxTxxx) = size(tmpdTxxx)
    (nzTxxz,nxTxxz) = size(tmpdTxxz)
    bxTxx = exp.(-(tmpdTxxx/kx.+ax)*medium.dt).-1
    bzTxx = exp.(-(tmpdTxxz/kz.+az)*medium.dt).-1
    axTxx = tmpdTxxx./(kx*(tmpdTxxx+kx*ax)).*bxTxx
    azTxx = tmpdTxxz./(kz*(tmpdTxxz+kz*az)).*bzTxx
    PTxxBVx = spzeros(nwf.BDntpp[1],nwf.BDntpp[2])
    PTxxBVz = spzeros(nwf.BDntpp[1],nwf.BDntpp[2])


    # Tzz
    (nz,nx) = nwf.ntpp
    tmpdTzzz,tmpdTzzx,ax,az = DampCoeff(nz, nx, medium.nDZ, medium.nHX, ext, maxv, iflag)
    (nzTzzx,nxTzzx) = size(tmpdTzzx)
    (nzTzzz,nxTzzz) = size(tmpdTzzz)
    bxTzz = exp.(-(tmpdTzzx/kx.+ax)*medium.dt).-1
    bzTzz = exp.(-(tmpdTzzz/kz.+az)*medium.dt).-1
    axTzz = tmpdTzzx./(kx*(tmpdTzzx+kx*ax)).*bxTzz
    azTzz = tmpdTzzz./(kz*(tmpdTzzz+kz*az)).*bzTzz
    PTzzBVx = spzeros(nwf.BDntpp[1],nwf.BDntpp[2])
    PTzzBVz = spzeros(nwf.BDntpp[1],nwf.BDntpp[2])

    # Txz
    (nz,nx) = nwf.ntxz
    tmpdTxzz,tmpdTxzx,ax,az = DampCoeff(nz, nx, medium.nDZ, medium.nHX, ext, maxv, iflag)
    (nzTxzx,nxTxzx) = size(tmpdTxzx)
    (nzTxzz,nxTxzz) = size(tmpdTxzz)
    bxTxz = exp.(-(tmpdTxzx/kx.+ax)*medium.dt).-1
    bzTxz = exp.(-(tmpdTxzz/kz.+az)*medium.dt).-1
    axTxz = tmpdTxzx./(kx.*(tmpdTxzx.+kx*ax)).*bxTxz
    azTxz = tmpdTxzz./(kz.*(tmpdTxzz.+kz*az)).*bzTxz
    PTxzBVx = spzeros(nwf.BDntxz[1],nwf.BDntxz[2])
    PTxzBVz = spzeros(nwf.BDntxz[1],nwf.BDntxz[2])

    return BD2d(bxVx, bzVx, axVx, azVx,
                      bxVz, bzVz, axVz, azVz,
                      bxTxx, bzTxx, axTxx, azTxx,
                      bxTzz, bzTzz, axTzz, azTzz,
                      bxTxz, bzTxz, axTxz, azTxz,
                      PVxBTxx, PVxBTxz,
                      PVzBTzz, PVzBTxz,
                      PTxxBVx, PTxxBVz,
                      PTzzBVx, PTzzBVz,
                      PTxzBVx, PTxzBVz)


end

# This function introduce Damping Coefficient to each damping layer
function DampCoeff(nz::Int64, nx::Int64, nDZ::Int64, nHX::Int64,
    ext::Int64, maxv::Float64, iflag::Int64)

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
        ax = ones(nz+2*ext,nx+2*ext)
        az = ones(nz+2*ext,nx+2*ext)

        if (nx==nHX)

            for i = 1 : ext
                dWFx[:,i] = -maxv/ext * log(R) * (m*(ext+1-i )/ext + n*((ext+1-i )/ext)^ll)
                ax[:,i] = ax[:,i].*tmp_a[i]
            end

            for i = nx+ext+1 : nx+2*ext
                dWFx[:,i] = -maxv/ext * log(R) * (m*(i-nx-ext)/ext + n*((i-nx-ext)/ext)^ll)
                ax[:,i] = ax[:,i].*tmp_a[nx+2*ext-i+1]
            end

        elseif (nx+1==nHX)

            for i = 1 : ext
                dWFx[:,i] = -maxv/ext * log(R) * (m*(ext+0.5-i)/ext + n*((ext+0.5-i)/ext)^ll)
                ax[:,i] = ax[:,i].*tmp_a[i]
            end
            for i = (nx+ext+1) : (nx+2*ext)
                dWFx[:,i] = -maxv/ext * log(R) * (m*(i-(nx)-ext-0.5)/ext + n*((i-(nx)-ext-0.5)/ext)^ll)
                ax[:,i] = ax[:,i].*tmp_a[nx+2*ext-i+1]
            end

        else error("please check your wavefield size")

        end

        if (nz==nDZ)

            for i = 1 : ext
                dWFz[i,:] = -maxv/ext * log(R) * (m*(ext+1-i )/ext + n*((ext+1-i )/ext)^ll)
                az[i,:] = az[i,:].*tmp_a[i]
            end
            for i = (nz+ext+1) : nz+2*ext
                dWFz[i,:] = -maxv/ext * log(R) * (m*(i-nz-ext)/ext + n*((i-nz-ext)/ext)^ll)
                az[i,:] = az[i,:].*tmp_a[nz+2*ext-i+1]
            end

        elseif (nz+1==nDZ)

            for i = 1 : ext
                dWFz[i,:] = -maxv/ext * log(R) * (m*(ext+0.5-i)/ext + n*((ext+0.5-i)/ext)^ll)
                az[i,:] = az[i,:].*tmp_a[i]
            end
            for i = nz+ext+1 : nz+2*ext
                dWFz[i,:] = -maxv/ext * log(R) * (m*(i-(nz)-ext-0.5)/ext + n*((i-(nz)-ext-0.5)/ext)^ll)
                az[i,:] = az[i,:].*tmp_a[nz+2*ext-i+1]
            end
            else error("please check your wavefield size")
        end

    elseif iflag == 1 # free surface

        dWFx = spzeros(nz+ext,nx+2*ext)
        dWFz = spzeros(nz+ext,nx+2*ext)
        ax = spzeros(nz+ext,nx+2*ext)
        az = spzeros(nz+ext,nx+2*ext)

        if (nx==nHX)

            for i = 1 : ext
                dWFx[:,i] = -maxv/ext * log(R) * (m*(ext+1-i )/ext + n*((ext+1-i )/ext)^ll)
                ax[:,i] = ax[:,i].*tmp_a[i]
            end

            for i = nx+ext+1 : nx+2*ext
                dWFx[:,i] = -maxv/ext * log(R) * (m*(i-nx-ext)/ext + n*((i-nx-ext)/ext)^ll)
                ax[:,i] = ax[:,i].*tmp_a[nx+2*ext-i+1]
            end

        elseif (nx+1==nHX)

            for i = 1 : ext
                dWFx[:,i] = -maxv/ext * log(R) * (m*(ext+0.5-i)/ext + n*((ext+0.5-i)/ext)^ll)
                ax[:,i] = ax[:,i].*tmp_a[i]
            end
            for i = (nx+ext+1) : (nx+2*ext)
                dWFx[:,i] = -maxv/ext * log(R) * (m*(i-(nx)-ext-0.5)/ext + n*((i-(nx)-ext-0.5)/ext)^ll)
                ax[:,i] = ax[:,i].*tmp_a[nx+2*ext-i+1]
            end
        else error("please check your wavefield size")
        end

        if (nz==nDZ)

            for i = (nz+1) : nz+ext
                dWFz[i,:] = -maxv/ext * log(R) * (m*(i-nz)/ext + n*((i-nz)/ext)^ll)
                az[i,:] = az[i,:].*tmp_a[nz+ext-i+1]
            end

        elseif (nz+1==nDZ)

            for i = nz+1 : nz+ext
                dWFz[i,:] = -maxv/ext * log(R) * (m*(i-(nz)-0.5)/ext + n*((i-(nz)-0.5)/ext)^ll)
                az[i,:] = az[i,:].*tmp_a[nz+ext-i+1]
            end
        else error("please check your wavefield size")
        end

    end

    return dWFz, dWFx, ax, az
end
