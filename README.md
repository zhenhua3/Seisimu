# Seisimu

This is a package for seismic wave propagation simulation in 2D/3D and Acoustic/Elastic medium using openmp parallel computing.

Package is not registered in Julia so to use the package one needs to add path of the package to Julia default searching path by adding the following line to ".juliarc.jl" in your home folder,
"Push!(LOAD_PATH,"where/you/store/the/package/src")".
Then when you can start use the package by typing "using seisimu" in julia pearl.


%%%%%%%%%%%%%%%%%%%%%%%%% Modelling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#=== Model initiation ===#

model = initmodel([Args]; [Optional Args])

Args in order: (1) P velocity, S velocity: 1D or text/bin file (m/s)

               (2) Density: Single value (kg/m3)

               (3) Depth, Horizon X: 2D medium / Depth, Horizon X, Horizon Y: 3D medium (m)

               (4) Peak Frequency: Single value (Hz)

               (5) Simulation time: Single value (s)

Optional Args: (6) PML boundary extension: [ext] = 5, 10, 20

               (7) Surface type: [iflag] = 1: free surface / 2: unlimited medium

               (8) Spatial grid intervel: [dz, dx]: 2D medium / [dz, dx, dy]: 3D medium (m)

               (9) Temporal sampling: [dt] (s)

               (10) Discrete simulation time: [nT]

Note: Words in square brackets of optional args are the key words when called, such "keyword = value".

#=== Source initiation ===#

source = initsource([Args])

Args in order: (1) Number of sources: Single value

               (2) 2D Array of the source locations: [z x]: 2D / [z x y]: 3D, coordinates are in columns

               (3) 1D Array of the source origin time: vector (s)

               (4) 2D Array of moment tensors: [m11 m22 m12]: 2D/ [m11 m22 m33 m12 m13 m23]: 3D

               (5) Waveforms of the sources: [vector(waveform1) vector(waveform2) ...]

               (6) Discrete model derived in model initiation

#=== Receiver initiation ===#

receiver = initrec([Args]])

Args in order: (1) Number of receivers: Single value

               (2) 2D Array of the receiver locations: [z x]: 2D / [z x y]: 3D, coordinates are in columns

               (3) Discrete model derived in model initiation

#=== Extrapolation ===#

runsimu([Args]; [Optional Args])

Args in order: (1) Discrete model derived in model initiation

               (2) Discrete source derived in source initiation

Optional Args: (1) Output recordings path: [RecDataPath] = String

               (2) Output wavefield path: [WFDataPath] = String

               (3) Output slices: can be any slice or a serial number, such as 1:20

               (4) Output component: [OptCpnt], choosen from ["P", "Vx", "Vz", "Txx", "Tzz", "Txz", "V", "Tii"]

Note: when either slices or recording path are not void, output component cannot be void.

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

>>>>>>> 7d1cb359e596d62551810e900a7e2582b1b687a1
