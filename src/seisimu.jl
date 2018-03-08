module seisimu

    # using PyCall
    # @pyimport matplotlib.pyplot as plt

    include("initiation/initiation.jl")
    include("extrap/extrap.jl")
    include("IO/IO.jl")
    include("initiation/ForCUDA/gpuchunk.jl")

    # functions and types

    export initmodel,Ricker,initsource,runsimu,
    gpuchunk,acmod2d,elmod3d,acmod3d,elmod2d,addsou!

end
