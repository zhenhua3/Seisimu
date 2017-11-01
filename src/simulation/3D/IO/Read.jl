#========================== Read recording or wavefield information ===============================#
type RecInfo
  RowNum :: Int64
  ColNum :: Int64
  Comp :: String # provide information of output data component
  RecNum :: Int64 # Number of Receivers
  RecLoc :: Array{Float64} # Receiver coordinates
  TimeSampNum :: Int64 # Total time sampling Tn
end

type WFInfo
  RowNum :: Int64
  ColNum :: Int64
  Comp :: String
  N_BDVx :: Array{Int64,1}
  N_BDVz :: Array{Int64,1}
  N_BDTpp :: Array{Int64,1}
  N_BDTxz :: Array{Int64,1}
  TimeSampNum :: Int64
  TimeSamp :: Array{Float64,1}
end

function ReadInfo(Path :: String, datatype :: String)
  fid = open(Path,"r")
  if datatype == "rec"
    iflag = read(fid,Int64,1)
    if iflag[1] == 1801
      comp = "[Vx; Vz]"
      compnum = 2
    elseif iflag[1] == 1802
      comp = "[P]"
      compnum = 1
    elseif iflag[1] == 1803
      comp = "[Txz]"
      compnum = 1
    elseif iflag[1] == 1807
      comp = "[Vx; Vz; P; Txz]"
      compnum = 4
    elseif iflag[1] == 1804
      comp = "[Vx; Vz; P]"
      compnum = 3
    elseif iflag[1] == 1805
      comp = "[Vx; Vz; Txz]"
      compnum = 3
    elseif iflag[1] == 1806
      comp = "[P; Txz]"
      compnum = 2
    else error("Recording data types are not defined")
    end
    recnum,timesampnum = read(fid,Int64,2) # Number of Receivers and Total time sampling Tn
    RowNum = compnum*recnum
    ColNum = timesampnum
    RecLoc = read(fid, Float64, 2*recnum)
    RecLoc = reshape(RecLoc, 2, recnum)
    close(fid)
    Info = RecInfo(RowNum, ColNum, comp,recnum,RecLoc,timesampnum)
    println(fieldnames(Info))
    return Info

  elseif datatype == "wf"
    iflag = read(fid,Int64,1)
    tppsize = read(fid,Int64,2)
    vxsize = [tppsize[1],tppsize[2]-1]
    vzsize = [tppsize[1]-1,tppsize[2]]
    txzsize = [tppsize[1]-1,tppsize[2]-1]
    timesampnum = read(fid,Int64,1)
    timesamp = read(fid,Float64,timesampnum[1])
    ColNum = timesampnum[1]

    if iflag[1] == 1801
      comp = "[Vx; Vz]"
      RowNum = vxsize[1]*vxsize[2]+vzsize[1]*vzsize[2]
    elseif iflag[1] == 1802
      comp = "[P]"
      RowNum = tppsize[1]*tppsize[2]
    elseif iflag[1] == 1803
      comp = "[Txz]"
      RowNum = txzsize[1]*txzsize[2]
    elseif iflag[1] == 1807
      comp = "[Vx; Vz; P; Txz]"
      RowNum = tppsize[1]*tppsize[2]+vxsize[1]*vxsize[2]+vzsize[1]*vzsize[2]+txzsize[1]*txzsize[2]
    elseif iflag[1] == 1804
      comp = "[Vx; Vz; P]"
      RowNum = tppsize[1]*tppsize[2]+vxsize[1]*vxsize[2]+vzsize[1]*vzsize[2]
    elseif iflag[1] == 1805
      comp = "[Vx; Vz; Txz]"
      RowNum = vxsize[1]*vxsize[2]+vzsize[1]*vzsize[2]+txzsize[1]*txzsize[2]
    elseif iflag[1] == 1806
      comp = "[P; Txz]"
      RowNum = tppsize[1]*tppsize[2]+txzsize[1]*txzsize[2]
    elseif iflag[1] == 1808
      comp = "[Vxx; Vxz; Vzx; Vzz; Txxx; Txxz; Tzzx; Tzzz; Txzx; Txzz]"
      RowNum = 2*(2*tppsize[1]*tppsize[2]+vxsize[1]*vxsize[2]+vzsize[1]*vzsize[2]+txzsize[1]*txzsize[2])
    else error("Recording data types are not defined")
    end
    close(fid)
    Info = WFInfo(RowNum,ColNum,comp,vxsize,vzsize,tppsize,txzsize,timesampnum[1],timesamp)
    print(fieldnames(Info))
    return Info

  else error("choose datatype from 'rec' (recordings) or 'wf' (wavefield) ")
  end
end

#========================== Read Data ===============================#
function ReadData(fid::IOStream, datatype::DataType, datalength::Int64)
  # This function read Recordings
  data = read(fid,datatype,datalength)
  rintln(eof(fid))
  return data
end

function ReadData(path::String, datatype::DataType, datalength::Int64)
  fid = open(path,"r")
  data = read(fid,datatype,datalength)
  println(eof(fid))
  return data
  close(fid)
end

function ReadData(path::String, datatype::DataType, RowNum::Int64, ColNum::Int64)
  fid = open(path,"r")
  data = read(fid,datatype,RowNum*ColNum)
  println(eof(fid))
  close(fid)
  data = reshape(data, RowNum, ColNum)
  return data
end

function ReadData(path::String, datatype::DataType)
  fid = open(path,"r")
  data = []
  while !eof(fid)
    tmp = read(fid,datatype,1)
    append!(data,tmp)
  end
  close(fid)
end

#========================== Read Simulation Parameter ===============================#
function ReadSimuPara(Path::String; SouFlag = false, mor = false, FDorder = 4)
  fid = open(Path,"r")
  iflag, ext = read(fid,Int64,2)
  modelsize = read(fid,Int64,2) # Model size Nz and Nx with PML boundary

  if iflag == 2 # unlimited medium
    NDep = modelsize[1] - 2*ext
    NHor = modelsize[2] - 2*ext
  elseif iflag == 1 # free surface
    NDep = modelsize[1] - ext
    NHor = modelsize[2] - 2*ext
  end

  PML_NDep = modelsize[1]
  PML_NHor = modelsize[2]
  PMLVP = read(fid,Float64,modelsize[1]*modelsize[2])
  PMLVP = reshape(PMLVP,modelsize[1],modelsize[2])
  PMLVS = read(fid,Float64,modelsize[1]*modelsize[2])
  PMLVS = reshape(PMLVS,modelsize[1],modelsize[2])
  PMLRho = read(fid,Float64,modelsize[1]*modelsize[2])
  PMLRho = reshape(PMLRho,modelsize[1],modelsize[2])
  Lambda = PMLVP.^2.*PMLRho - 2*PMLVS.^2.*PMLRho
  Mu = PMLVS.^2.*PMLRho
  dx,dz,dt = read(fid,Float64,3)
  Dep = NDep*dz
  Hor = NHor*dx
  Tn,sn = read(fid,Int64,2)
  TST = Tn*dt
  Medium = Model(PMLVP, PMLVS, PMLRho, Lambda, Mu, dx, dz, dt, Dep, Hor, NDep, NHor, PML_NDep, PML_NHor,TST,Tn)
  PF = read(fid,Float64,1)

  N_BDTpp = modelsize
  N_BDTxz = [N_BDTpp[1]-1; N_BDTpp[2]-1]
  N_BDVx = [N_BDTpp[1]; N_BDTpp[2]-1]
  N_BDVz = [N_BDTpp[1]-1; N_BDTpp[2]]
  if iflag == 1
    N_Tpp = [N_BDTpp[1]-ext; N_BDTpp[2]-2*ext]
    N_Txz = [N_BDTxz[1]-ext; N_BDTxz[2]-2*ext]
    N_Vx = [N_BDVx[1]-ext; N_BDVx[2]-2*ext]
    N_Vz = [N_BDVz[1]-ext; N_BDVz[2]-2*ext]
  elseif iflag == 2
    N_Tpp = [N_BDTpp[1]-2*ext; N_BDTpp[2]-2*ext]
    N_Txz = [N_BDTxz[1]-2*ext; N_BDTxz[2]-2*ext]
    N_Vx = [N_BDVx[1]-2*ext; N_BDVx[2]-2*ext]
    N_Vz = [N_BDVz[1]-2*ext; N_BDVz[2]-2*ext]
  end
  wfsize = WF2DSize(N_BDTpp,N_BDTxz,N_BDVx,N_BDVz,N_Tpp,N_Txz,N_Vx,N_Vz)

  if mor == true
    FDC = FDCoeff(FDorder)
    FD = Fdmtx(wfsize,Medium,FDC,ext,iflag)
  elseif mor == false
    FDC = FDCoeff(FDorder)
    FD = Fdmtx(wfsize,Medium,FDC,ext)
  end

  if SouFlag == true
      posn = read(fid,Int64,sn*2)
      posn = reshape(posn,2,sn)
      position_z = posn[1,:]*dz
      position_x = posn[2,:]*dx
      position = [position_z'; position_x']
      if sn == 1
        position = vec(position)
      end
      Waveform = read(fid,Float64,Tn*sn)
      Waveform = reshape(Waveform,Tn,sn)
      if sn == 1
        Waveform = vec(Waveform)
      end
      sotn = read(fid,Int64,sn)
      sot = sotn.*dt
      if sn == 1
        sot = sot[1]
      end
      stpcode = read(fid,Int64,sn)
      sourcetype = Array{String}(sn)
        for i = 1:sn
          if stpcode[i] == 1
            sourcetype[i] = "expl"
          elseif stpcode[i] == 2
            sourcetype[i] = "sfx"
          elseif stpcode[i] == 3
            sourcetype[i] = "sfz"
          elseif stpcode[i] == 4
            sourcetype[i] = "scx"
          elseif stpcode[i] == 5
            sourcetype[i] = "scz"
          elseif stpcode[i] == 6
            sourcetype[i] = "DC"
          end
        end
        if sn == 1
          sourcetype = sourcetype[1]
        end
      Sou = InitSou(sn,position,sot,sourcetype,Waveform,Medium,wfsize; ext = ext, iflag = iflag)
      Sou.waveform = Waveform
      close(fid)
      return Medium,Sou,PF,wfsize,FD
  elseif SouFlag == false
    close(fid)
    return Medium,PF,wfsize,FD
  end
end

#==================== Read Basis ===========================#
function ReadBasis(Path::String; ColNum = nothing)
  fid = open(Path,"r")
  if ColNum == nothing
    RowNumber, ColumNumber = read(fid,Int64,2)
    Q = read(fid,Float64,RowNumber*ColumNumber)
    Q = reshape(Q, RowNumber, ColumNumber)
    close(fid)
  elseif ColNum != nothing
    RowNumber, ColumNumber = read(fid,Int64,2)
    if ColNum == 0
      error("Column Number has to be larger than 0")
    elseif ColNum > ColumNumber
      error("Column Nuber has to be equal to or smaller than $ColumNumber")
    else
      Q = read(fid,Float64,RowNumber*ColNum)
      Q = reshape(Q, RowNumber, ColNum)
      close(fid)
    end
  end
  return Q
end
