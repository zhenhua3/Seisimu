type dpcf
    bfull::Array{Float32,1}
    bhalf::Array{Float32,1}
    afull::Array{Float32,1}
    ahalf::Array{Float32,1}
end


function DampCoeff(medium::Union{elastic2d,acoustic2d,elastic3d,acoustic3d})
    ext = medium.ext
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
    tmp_a = linspace(a_max,0,ext)
    dWFfull = zeros(ext)
    dWFhalf = zeros(ext)
    maxv = maximum(medium.pvel)
    for i = 1 : ext
        dWFfull[i] = -maxv/ext * log(R) * (m*(ext+1-i)/ext + n*((ext+1-i )/ext)^ll)
        dWFhalf[i] = -maxv/ext * log(R) * (m*(ext+0.5-i)/ext + n*((ext+0.5-i)/ext)^ll)
    end
    bhalf = Float32.(exp.(-(dWFhalf.+tmp_a).*medium.dt))
    bfull = Float32.(exp.(-(dWFfull.+tmp_a).*medium.dt))
    ahalf = Float32.(dWFhalf./(dWFhalf.+tmp_a).*(bhalf.-1))
    afull = Float32.(dWFfull./(dWFfull.+tmp_a).*(bfull.-1))
    return dpcf(bfull,bhalf,afull,ahalf)
end


            #==========================#
            #=== Obsorbing Boundary ===#
            #==========================#
#=== elastic 2d ===#
function BD(nwf::nelwf2d, medium::elastic2d)

    ext = medium.ext
    iflag = medium.iflag
    dpcf = DampCoeff(medium)

    # Vx
    indexZ = nwf.BDnvx[1]
    indexX = nwf.BDnvx[2] - nwf.nvx[2]
    PVxBTxx = zeros(Float64,indexZ,indexX)

    indexZ = nwf.BDnvx[1] - nwf.nvx[1]
    indexX = nwf.BDnvx[2]
    PVxBTxz = zeros(Float64,indexZ,indexX)

    # Vz
    indexZ = nwf.BDnvz[1] - nwf.nvz[1]
    indexX = nwf.BDnvz[2]
    PVzBTzz = zeros(Float64,indexZ,indexX)

    indexZ = nwf.BDnvz[1]
    indexX = nwf.BDnvz[2] - nwf.nvz[2]
    PVzBTxz = zeros(Float64,indexZ,indexX)

    # Tzz & Txx
    indexZ = nwf.BDntpp[1] - nwf.ntpp[1]
    indexX = nwf.BDntpp[2]
    PTzzBVz = zeros(Float64,indexZ,indexX)
    PTxxBVz = zeros(Float64,indexZ,indexX)

    indexZ = nwf.BDntpp[1]
    indexX = nwf.BDntpp[2] - nwf.ntpp[2]
    PTzzBVx = zeros(Float64,indexZ,indexX)
    PTxxBVx = zeros(Float64,indexZ,indexX)

    # Txz
    indexZ = nwf.BDntxz[1] - nwf.ntxz[1]
    indexX = nwf.BDntxz[2]
    PTxzBVx = zeros(Float64,indexZ,indexX)

    indexZ = nwf.BDntxz[1]
    indexX = nwf.BDntxz[2] - nwf.ntxz[2]
    PTxzBVz = zeros(Float64,indexZ,indexX)

    return elbd2d(Float64.(dpcf.bfull), Float64.(dpcf.bhalf), Float64.(dpcf.afull), Float64.(dpcf.ahalf),
                      PVxBTxx, PVxBTxz,
                      PVzBTzz, PVzBTxz,
                      PTxxBVx, PTxxBVz,
                      PTzzBVx, PTzzBVz,
                      PTxzBVx, PTxzBVz)
end


#=== acoustic 2d ===#
function BD(nwf::nacwf2d, medium::acoustic2d)
    ext = medium.ext
    iflag = medium.iflag
    dpcf = DampCoeff(medium)

    # Vx
    indexZ = nwf.BDnvx[1]
    indexX = nwf.BDnvx[2] - nwf.nvx[2]
    PVxBTxx = zeros(Float64,indexZ,indexX)

    # Vz
    indexZ = nwf.BDnvz[1] - nwf.nvz[1]
    indexX = nwf.BDnvz[2]
    PVzBTzz = zeros(Float64,indexZ,indexX)

    # Tzz & Txx
    indexZ = nwf.BDntpp[1] - nwf.ntpp[1]
    indexX = nwf.BDntpp[2]
    PTzzBVz = zeros(Float64,indexZ,indexX)
    PTxxBVz = zeros(Float64,indexZ,indexX)

    indexZ = nwf.BDntpp[1]
    indexX = nwf.BDntpp[2] - nwf.ntpp[2]
    PTzzBVx = zeros(Float64,indexZ,indexX)
    PTxxBVx = zeros(Float64,indexZ,indexX)

    return acbd2d(Float64.(dpcf.bfull), Float64.(dpcf.bhalf), Float64.(dpcf.afull), Float64.(dpcf.ahalf),
                      PVxBTxx, PVzBTzz, 
                      PTxxBVx, PTxxBVz,
                      PTzzBVx, PTzzBVz)
end





#=== elastic 3d ===#
function BD(nwf::nelwf3d, medium::elastic3d)

        ext = medium.ext
        iflag = medium.iflag
        dpcf = DampCoeff(medium)

        # Vx
        indexZ = nwf.BDnvx[1]
        indexX = nwf.BDnvx[2] - nwf.BDnvx[2]
        indexY = nwf.BDnvx[3]
        PVxBTxx = zeros(Float32,indexZ,indexX,indexY)

        indexZ = nwf.BDnvx[1] - nwf.nvx[1]
        indexX = nwf.BDnvx[2]
        indexY = nwf.BDnvx[3]
        PVxBTxz = zeros(Float32,indexZ,indexX,indexY)

        indexZ = nwf.BDnvx[1]
        indexX = nwf.BDnvx[2]
        indexY = nwf.BDnvx[3] - nwf.nvx[3]
        PVxBTxy = zeros(Float32,indexZ,indexX,indexY)

        # Vy
        indexZ = nwf.BDnvy[1]
        indexX = nwf.BDnvy[2] - nwf.BDnvx[2]
        indexY = nwf.BDnvy[3]
        PVyBTxy = zeros(Float32,indexZ,indexX,indexY)

        indexZ = nwf.BDnvy[1] - nwf.nvy[1]
        indexX = nwf.BDnvy[2]
        indexY = nwf.BDnvy[3]
        PVyBTyz = zeros(Float32,indexZ,indexX,indexY)

        indexZ = nwf.BDnvx[1]
        indexX = nwf.BDnvx[2]
        indexY = nwf.BDnvx[3] - nwf.nvx[3]
        PVyBTyy = zeros(Float32,indexZ,indexX,indexY)

        # Vz
        indexZ = nwf.BDnvz[1] - nwf.nvz[1]
        indexX = nwf.BDnvz[2]
        indexY = nwf.BDnvz[3]
        PVzBTzz = zeros(Float32,indexZ,indexX,indexY)

        indexZ = nwf.BDnvz[1]
        indexX = nwf.BDnvz[2] - nwf.nvz[2]
        indexY = nwf.BDnvz[3]
        PVzBTxz = zeros(Float32,indexZ,indexX,indexY)

        indexZ = nwf.BDnvz[1]
        indexX = nwf.BDnvz[2]
        indexY = nwf.BDnvz[3] - nwf.nvz[3]
        PVzBTyz = zeros(Float32,indexZ,indexX,indexY)


        # Tzz & Txx & Tyy
        indexZ = nwf.BDntpp[1] - nwf.ntpp[1]
        indexX = nwf.BDntpp[2]
        indexY = nwf.BDntpp[3]
        PTzzBVz = zeros(Float32,indexZ,indexX,indexY)
        PTxxBVz = zeros(Float32,indexZ,indexX,indexY)
        PTyyBVz = zeros(Float32,indexZ,indexX,indexY)

        indexZ = nwf.BDntpp[1]
        indexX = nwf.BDntpp[2] - nwf.ntpp[2]
        indexY = nwf.BDntpp[3]
        PTzzBVx = zeros(Float32,indexZ,indexX,indexY)
        PTxxBVx = zeros(Float32,indexZ,indexX,indexY)
        PTyyBVx = zeros(Float32,indexZ,indexX,indexY)

        indexZ = nwf.BDntpp[1]
        indexX = nwf.BDntpp[2]
        indexY = nwf.BDntpp[3] - nwf.ntpp[3]
        PTzzBVy = zeros(Float32,indexZ,indexX,indexY)
        PTxxBVy = zeros(Float32,indexZ,indexX,indexY)
        PTyyBVy = zeros(Float32,indexZ,indexX,indexY)

        # Txz
        indexZ = nwf.BDntxz[1] - nwf.ntxz[1]
        indexX = nwf.BDntxz[2]
        indexY = nwf.BDntxz[3]
        PTxzBVx = zeros(Float32,indexZ,indexX,indexY)

        indexZ = nwf.BDntxz[1]
        indexX = nwf.BDntxz[2] - nwf.ntxz[2]
        indexY = nwf.BDntxz[3]
        PTxzBVz = zeros(Float32,indexZ,indexX,indexY)

        # Tyz
        indexZ = nwf.BDntyz[1] - nwf.ntyz[1]
        indexX = nwf.BDntyz[2]
        indexY = nwf.BDntyz[3]
        PTyzBVy = zeros(Float32,indexZ,indexX,indexY)

        indexZ = nwf.BDntyz[1]
        indexX = nwf.BDntyz[2]
        indexY = nwf.BDntyz[3] - nwf.ntyz[3]
        PTyzBVz = zeros(Float32,indexZ,indexX,indexY)

        # Txy
        indexZ = nwf.BDntxy[1]
        indexX = nwf.BDntxy[2] - nwf.ntxy[2]
        indexY = nwf.BDntxy[3]
        PTxyBVy = zeros(Float32,indexZ,indexX,indexY)

        indexZ = nwf.BDntxy[1]
        indexX = nwf.BDntxy[2]
        indexY = nwf.BDntxy[3] - nwf.ntxy[3]
        PTxyBVx = zeros(Float32,indexZ,indexX,indexY)

    return elbd3d(dpcf.bfull,dpcf.bhalf,dpcf.afull,dpcf.ahalf,
                  PVzBTzz, PVzBTxz, PVzBTyz,
                  PVxBTxx, PVxBTxz, PVxBTxy,
                  PVyBTyy, PVyBTyz, PVyBTxy,
                  PTzzBVz, PTzzBVx, PTzzBVy,
                  PTxxBVz, PTxxBVx, PTxxBVy,
                  PTyyBVz, PTyyBVx, PTyyBVy,
                  PTxzBVx, PTxzBVz,
                  PTyzBVz, PTyzBVy,
                  PTxyBVx, PTxyBVy)
end


#=== acoustic ===#
function BD(nwf::nacwf3d, medium::acoustic3d)

    ext = medium.ext
    iflag = medium.iflag
    dpcf = DampCoeff(medium)

    # Vx
    indexZ = nwf.BDnvx[1]
    indexX = nwf.BDnvx[2] - nwf.BDnvx[2]
    indexY = nwf.BDnvx[3]
    PVxBTxx = zeros(Float32,indexZ,indexX,indexY)

    # Vy
    indexZ = nwf.BDnvx[1]
    indexX = nwf.BDnvx[2]
    indexY = nwf.BDnvx[3] - nwf.nvx[3]
    PVyBTyy = zeros(Float32,indexZ,indexX,indexY)

    # Vz
    indexZ = nwf.BDnvz[1] - nwf.nvz[1]
    indexX = nwf.BDnvz[2]
    indexY = nwf.BDnvz[3]
    PVzBTzz = zeros(Float32,indexZ,indexX,indexY)

    # Tzz & Txx & Tyy
    indexZ = nwf.BDntpp[1] - nwf.ntpp[1]
    indexX = nwf.BDntpp[2]
    indexY = nwf.BDntpp[3]
    PTzzBVz = zeros(Float32,indexZ,indexX,indexY)
    PTxxBVz = zeros(Float32,indexZ,indexX,indexY)
    PTyyBVz = zeros(Float32,indexZ,indexX,indexY)

    indexZ = nwf.BDntpp[1]
    indexX = nwf.BDntpp[2] - nwf.ntpp[2]
    indexY = nwf.BDntpp[3]
    PTzzBVx = zeros(Float32,indexZ,indexX,indexY)
    PTxxBVx = zeros(Float32,indexZ,indexX,indexY)
    PTyyBVx = zeros(Float32,indexZ,indexX,indexY)

    indexZ = nwf.BDntpp[1]
    indexX = nwf.BDntpp[2]
    indexY = nwf.BDntpp[3] - nwf.ntpp[3]
    PTzzBVy = zeros(Float32,indexZ,indexX,indexY)
    PTxxBVy = zeros(Float32,indexZ,indexX,indexY)
    PTyyBVy = zeros(Float32,indexZ,indexX,indexY)

    return acbd3d(dpcf.bfull,dpcf.bhalf,dpcf.afull,dpcf.ahalf,
                  PVzBTzz, PVxBTxx, PVyBTyy,
                  PTzzBVz, PTzzBVx, PTzzBVy,
                  PTxxBVz, PTxxBVx, PTxxBVy,
                  PTyyBVz, PTyyBVx, PTyyBVy)
end
