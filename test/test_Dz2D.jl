function test_Dz2D()

    include("/home/lzh/Project/STGDFD/src/FDMtx/FDCoeff.jl")
    include("/home/lzh/Project/STGDFD/src/FDMtx/GFDMtx.jl")
    include("/home/lzh/Project/STGDFD/src/PDiff/2D/Dz2D.jl")

    path = "~/Desktop/PROPAGATOR/test.bin"

    fid = open("/home/lzh/Desktop/PROPAGATOR/test.bin","r")
    test = read(fid, Float32, 1600)
    close(fid)

    test = reshape(test,(40,40))

    fid = open("/home/lzh/Desktop/PROPAGATOR/test_z.bin","r")
    test_z = read(fid, Float32, 1600)
    close(fid)

    test_z = reshape(test_z,(40,40))

    Nz = 41

    FDC = FDCoeff(4)

    GFDM = GFDMtx(40,FDC)

    (DzWF, FDM) = Dz2D(test, Nz, GFDM, FDC)

    norm(DzWF[3:38,3:38] - test_z[3:38,3:38])
end





