
    using seisimu

    pvel = 3500;
    rho = 2.8;
    DZ = 4000;
    HX = 4000;
    HY = 4000;
    pkf = 50;
    T = 2;
    ext = 10;
    iflag = 2;
    dz = nothing;
    dx = nothing;
    dt= nothing;
    nT= nothing;
    mor = false;

    model = init3DAC(pvel,rho,DZ,HX,HY,pkf,T);
    path = "example_gpu.bin"
    intvl = 3
    threadim = [16,8,8]

    # acgroup = MemCpyGroups(model,10)
    # # #
    # acgroup = acMemcpy(model,10)

# @time run!(model, path; intvl=20, threadim = threadim, AssignedStreamNum = 6)

@time begin
    for i in 1:10#model.medium.nT
        run!(model)
    end
end
