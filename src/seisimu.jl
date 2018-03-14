module seisimu

    # using PyCall
    # @pyimport matplotlib.pyplot as plt

    include("initiation/initiation.jl")
    include("extrap/run_base.jl")
    include("extrap/run_simulation.jl")
    include("IO/IO.jl")

    # GPU part is currently unenabled
    # include("initiation/ForCUDA/gpuchunk.jl")

    # functions and types
    export initmodel,Ricker,initsource,runsimu,
    acmod2d,elmod3d,acmod3d,elmod2d,addsou!,run!

end
