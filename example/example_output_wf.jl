# function Modelling example

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# In this example, a simple homogeneous elastic modelling is shown with proper#
# settings. The example is a multisource example.                             #
# To run the example, copy the following code into Julia interface. You may   #
# need to check if the path has been added
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

    using seisimu

    pvel = 3500;
    svel = nothing;
    rho = 2.8;
    DZ = 1000;
    HX = 1000;
    pkf = 50;
    T = 0.5;
    ext = 10;
    iflag = 2;
    dz = nothing;
    dx = nothing;
    dt= nothing;
    nT= nothing;
    mor = false;

    model = initmodel(pvel,svel,rho,DZ,HX,pkf,T);

    sn = 2;
    loc = [500 500; 800 300];
    ot = [0.05, 0.08];
    tp = ["expl", "expl"];
    wavelet = Ricker(pkf,model.medium.dt);

    waveform = [wavelet wavelet];

    sou = initsource(sn,loc,ot,tp,waveform,model);

    OptCpnt = ["V"];
    slices = 50:50:1000;
    wfpath = "wavefields.bin";

    writeinfo(mode,sou,slices,OptCpnt)
    @time extrap2d(model,sou; Slices = slices, OptCpnt = OptCpnt, WFDataPath = wfpath)
