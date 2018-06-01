
    using seisimu

    pvel = 3500;
    rho = 2.8;
    DZ = 1000;
    HX = 1000;
    HY = 1000;
    pkf = 50;
    T = 0.4;
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

    run!(path, intvl, model, threadim)

@time begin
    for i in 1:model.medium.nT
        run!(model)
    end
end
