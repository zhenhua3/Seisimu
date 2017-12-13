function readmodel2d(ParaPath::String; mor = false)
    fid = open(ParaPath,"r")
    iflag, ext = read(fid,Int64,2)
    modelsize = read(fid,Int64,2) # Model size Nz and Nx with PML boundary

    if iflag == 2 # unlimited medium
        nDZ = modelsize[1] - 2*ext
        nHX = modelsize[2] - 2*ext
    elseif iflag == 1 # free surface
        nDZ = modelsize[1] - ext
        nHX = modelsize[2] - 2*ext
    end
    BDnDZ = modelsize[1]
    BDnHX = modelsize[2]

    pvel = read(fid,Float64,modelsize[1]*modelsize[2])
    pvel = reshape(pvel,modelsize[1],modelsize[2])
    svel = read(fid,Float64,modelsize[1]*modelsize[2])
    svel = reshape(svel,modelsize[1],modelsize[2])
    rho = read(fid,Float64,modelsize[1]*modelsize[2])
    rho = reshape(rho,modelsize[1],modelsize[2])
    lambda = pvel.^2.*rho - 2*svel.^2.*rho
    mu = svel.^2.*rho
    dx,dz,dt = read(fid,Float64,3)
    DZ = nDZ*dz
    HX = nHX*dx
    nT,sn = read(fid,Int64,2)
    T = nT*dt
    pkf = read(fid,Float64,1)
    medium = medium2d(pvel, svel, rho, lambda, mu, dx, dz, dt, DZ, HX,
        nDZ, nHX, BDnDZ, BDnHX,T,nT,pkf[1],ext,iflag)

    BDntpp = modelsize
    BDntxz = [BDntpp[1]-1, BDntpp[2]-1]
    BDnvx = [BDntpp[1], BDntpp[2]-1]
    BDnvz = [BDntpp[1]-1, BDntpp[2]]
    if iflag == 1
        ntpp = [BDntpp[1]-ext BDntpp[2]-2*ext]
        ntxz = [BDntxz[1]-ext, BDntxz[2]-2*ext]
        nvx = [BDnvx[1]-ext, BDnvx[2]-2*ext]
        nvz = [BDnvz[1]-ext, BDnvz[2]-2*ext]
    elseif iflag == 2
        ntpp = [BDntpp[1]-2*ext, BDntpp[2]-2*ext]
        ntxz = [BDntxz[1]-2*ext, BDntxz[2]-2*ext]
        nvx = [BDnvx[1]-2*ext, BDnvx[2]-2*ext]
        nvz = [BDnvz[1]-2*ext, BDnvz[2]-2*ext]
    end
    nwf = nwf2d(BDntpp,BDntxz,BDnvx,BDnvz,ntpp,ntxz,nvx,nvz)
    wf = initwf(nDZ,nHX,ext,iflag,mor)

    if mor == true
        FDC = FDCoeff(4)
        fd = fdmtx(nwf,medium,FDC,ext,iflag)
        return spmod2d(medium, wf, nwf, fd)
    elseif mor == false
        FDC = FDCoeff(4)
        fd = fdmtx(nwf,medium,FDC,ext)
        pml = BD(nwf, medium, ext, iflag)
        return nspmod2d(medium, wf, nwf, fd, pml)
    end
    close(fid)
end
