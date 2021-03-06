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
    DZ = 2000;
    HX = 2000;
    pkf = 30;
    T = 0.5;
    ext = 10;
    iflag = 2;
    dz = nothing;
    dx = nothing;
    dt= nothing;
    nT= nothing;
    mor = false;

    #=== elastic ===#
    model = initmodel(pvel,svel,rho,DZ,HX,pkf,T);

    #=== acoustic ===#
    model = initmodel(pvel,rho,DZ,HX,pkf,T);

    sn = 2;
    loc = [500 500; 800 300];
    ot = [0.05, 0.08];
    tp = ["expl", "expl"];
    wavelet = Ricker(pkf,model.medium.dt);

    waveform = [wavelet wavelet];

    sou = initsource(sn,loc,ot,tp,waveform,model);

    recloc = [1000 1000; 1100 1000; 1200 1000; 1300 1000];
    OptCpnt = ["Vx"];
    rec = initrec(4,recloc,OptCpnt,model);
    recpath = "recordings.bin";


    writeinfo(model,sou,rec,OptCpnt);
    @time extrap2d(model,sou; rec = rec, RecDataPath = recpath, OptCpnt = OptCpnt)
