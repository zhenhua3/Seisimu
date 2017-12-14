module seisimu

    # using PyCall
    # @pyimport matplotlib.pyplot as plt

    include("initiation/initiation.jl")
    include("extrap/extrap.jl")
    include("IO/IO.jl")

    # functions and types

    export initmodel,Ricker,initsource,runsimu

end
