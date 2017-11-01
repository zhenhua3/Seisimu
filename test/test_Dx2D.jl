function test_Dx2D()

    include("/home/lzh/Project/STGDFD/src/FiniteDiff/FDCoeff.jl")
    include("/home/lzh/Project/STGDFD/src/FiniteDiff/2D/Dx2D.jl")
    include("/home/lzh/Project/STGDFD/src/FiniteDiff/2D/Dz2D.jl")

    test = rand(11,12)

    Nx = 12;
    Nz = 11;

    FDC = FDCoeff(4);

    Dz= Dz2D(test, Nz, Nx, FDC);

end





