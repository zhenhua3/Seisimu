module SeiSimu

    include("simulation/2D/2D.jl")
    include("IO/IO.jl")

    # functions and types

    export extrap2d, Ricker,
    initsource, initrec,
    initmodel, writeinfo,
    readmodel2d, readrec, readwfinfo,
    readwfdata, nspmod2d, spfdmtx2d, nspfdmtx2d

end
