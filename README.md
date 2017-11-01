# SeismicBox
This is a package for seismic wave propagation simulation, migration and inversion for 2D/3D and Acoustic/Elastic medium

%%%%%%%%%%%%%%%%%%%%%%%%% Modelling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Staggered grid finite difference with unsplit PML boundary.
Functions that needed in the 2D simulation are introduced in the following content.

#=== Model initiation ===#
model = initmodel(v_p, v_s, rho, Depth, Horizon, PeakFreq, T;
                ext = 10, iflag = 2,
                dz = nothing, dx = nothing, dt= nothing,
                Tn= nothing, mor = false)

#=== Source initiation ===#
source = initsource(sn, position, origintime, sourcetype, waveform, model)

Source types "expl": Explosive source;
              "sfx": Single force in X direction;
              "sfz": Single force in Z direction;
              "scx": Single couple in X direction;
              "scz": Single couple in Z direction;
               "DC": Double couple source.

#=== Receiver initiation ===#
receiver = initrec(rn,receivers,OptCpnt,model)

#=== 2D Extrapolation ===#
Extrap2D(model,sou; rec = nothing, RecDataPath = nothing,
             WFDataPath = nothing, Slices = nothing,
                OptCpnt = nothing)

rec ::receiver
RecDataPath:: String
WFDataPath:: String
Slices::

#=== Data I/O ===#
writeinfo()

  ReadInfo() - Function can read wavefields or recordings information

    (1) Info = ReadInfo(Path, datatype)

    Note: 1. datatype includes two choices "rec" and "wf". "rec" is recordings and "wf" is wavefields.
            Target file is "info.bin"

  ReadSimuPara() - Function can read simulation parameter in bin file, then provide necessary data for remodelling

    (1) Medium, Sou, PeakFreq, wfsize, FD = ReadSimuPara(Path; SouFlag = true)
    (2) Medium, PeakFreq, wfsize, FD = ReadSimuPara(Path)

    Note: 1. parameter "SouFlag" has value: true or false. If true, output: Medium::Model, Sou::Union{Source, Array{Source}}, PeakFreq::Float64, wfsize::WF2DSize. If false, output: Medium::Model, PeakFreq::Float64, wfsize::WF2DSize.

  ReadData() - Function
    (1) Read data only with Input IOStream and known length
      data = ReadData(fid::IOStream, datatype, datalength)

    (2) Read data only with Input file path and known length
      data = ReadData(Path, datatype, datalength)

    (3) Read data and reshape data to a matrix and known length
      data = ReadData(Path, datatype, RowNum, ColNum)

    (4) Read data only but data length is unknown
      data = ReadData(Path, datatype)

      Note: datatype: Float64, Int64, etc
  ReadBasis() - Function

    (1) Read basis file and reshape basis to a certain size
      basis = ReadBasis(Path)
