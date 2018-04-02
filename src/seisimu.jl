module seisimu

    include("initiation/initiation.jl")
    include("extrap/run_base.jl")
    include("extrap/run_simulation.jl")
    include("IO/IO.jl")

    # GPU part is currently disabled
    # include("initiation/ForCUDA/gpuchunk.jl")

    # functions and types
    export init2DEL, init2DAC, init3DEL, init3DAC, reset!,
    Ricker,initMTsource,initSFsource,
    runsimu,
    acmod2d,elmod3d,acmod3d,elmod2d,
    addmt!,addsf!,run!

end
