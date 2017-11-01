# All Output parameters are turned into a number code in this function as follows
# Opt["V"] == true && Opt["P"] == false && Opt["Txz"] == false : 1801
# Opt["V"] == false && Opt["P"] == true && Opt["Txz"] == false : 1802
# Opt["V"] == false && Opt["P"] == false && Opt["Txz"] == true : 180
# Opt["V"] == true && Opt["P"] == true && Opt["Txz"] == false : 1804
# Opt["V"] == true && Opt["P"] == false && Opt["Txz"] == true : 1805
# Opt["V"] == false && Opt["P"] == true && Opt["Txz"] == true : 1806
# Opt["V"] == true && Opt["P"] == true && Opt["Txz"] == true : 1807
# Number code are embeded in Output data so that through reading function no need to
# specify input data type
# Opt is a short for Output
#======================= Simulate Multiple Source Seismic Recordings ==============================#
function Simulation{T,T1<:Real}(PF::T, WF1::WF2D,WF2::WF2D,WFSize::WF2DSize, Medium::Model, source::Array{Source},
  FD::NSFDMtx, pml::dampCoef, sn::Int64, Opt, receivers::Array{T1},rn::Int64; ext = 10, iflag = 2,
  infotxtpath = nothing, infobinpath = nothing, path_rec = "rec.bin", path_recinfo = "recinfo.bin")

  receivers = Array{Float64}(receivers)
  Rec, PML_Rec = InitRec(rn,receivers,Medium.dz,Medium.dx, ext, iflag)
  Tn = Medium.Tn
  Writeinfo(WFSize,source,sn,PF,iflag,Medium,PML_Rec,infotxtpath,infobinpath,Opt,ext)

  if path_rec[end-2:end] != "bin" || path_recinfo[end-2:end] != "bin"
    error("Output recording data has to be in 'bin' format")
  end

  if Opt["V"] == true && Opt["P"] == false && Opt["Txz"] == false

    #========= write Recordings info =========#
    fid_recinfo = open(path_recinfo, "a+")
    write(fid_recinfo,1801)
    write(fid_recinfo,rn) # Total Receiver Number
    write(fid_recinfo,Medium.Tn) # Total Time Sampling Tn
    write(fid_recinfo,receivers)
    close(fid_recinfo)
    #========== write Recordings =========#
    RecLoc_Vx = (PML_Rec[2,:]-1).*WFSize.N_BDVx[1] + PML_Rec[1,:]
    RecLoc_Vz = (PML_Rec[2,:]-1).*WFSize.N_BDVz[1] + PML_Rec[1,:]
    fid_record = open(path_rec, "a+")
    for it = 1 : Tn
      Addsource!(WF1, sn, source, it)
      Onestep2D!(WF1, WF2, FD, pml)
      OptRec = [WF1.VecBDVx[RecLoc_Vx];WF1.VecBDVz[RecLoc_Vz]]
      write(fid_record,OptRec)
    end
    close(fid_record)

  elseif Opt["V"] == false && Opt["P"] == true && Opt["Txz"] == false

    #========= write Recordings info =========#
    fid_recinfo = open(path_recinfo, "a+")
    write(fid_recinfo,1802)
    write(fid_recinfo,rn) # Total Receiver Number
    write(fid_recinfo,Medium.Tn) # Total Time Sampling Tn
    write(fid_recinfo,receivers)
    close(fid_recinfo)
    #=========== write Recordings ===========#
    RecLoc_P = (PML_Rec[2,:]-1).*WFSize.N_BDTpp[1] + PML_Rec[1,:]
    fid_record = open(path_rec, "a+")
    for it = 1 : Tn
      Addsource!(WF1, sn, source, it)
      Onestep2D!(WF1, WF2, FD, pml)
      OptRec = (WF1.VecBDTxx[RecLoc_P]+WF1.VecBDTzz[RecLoc_P])/2
      write(fid_record,OptRec)
    end
    close(fid_record)

  elseif Opt["V"] == false && Opt["P"] == false && Opt["Txz"] == true

    #========= write Recordings info =========#
    fid_recinfo = open(path_recinfo, "a+")
    write(fid_recinfo,1803)
    write(fid_recinfo,rn) # Total Receiver Number
    write(fid_recinfo,Medium.Tn) # Total Time Sampling Tn
    write(fid_recinfo,receivers)
    close(fid_recinfo)
    #=========== write Recordings ===========#
    RecLoc_Txz = (PML_Rec[2,:]-1).*WFSize.N_BDTxz[1] + PML_Rec[1,:]
    fid_record = open(path_rec, "a+")
    for it = 1 : Tn
      Addsource!(WF1, sn, source, it)
      Onestep2D!(WF1, WF2, FD, pml)
      OptRec = WF1.VecBDTxz[RecLoc_Txz]
      write(fid_record,OptRec)
    end
    close(fid_record)

  elseif Opt["V"] == true && Opt["P"] == true && Opt["Txz"] == true

    #========= write Recordings info =========#
    fid_recinfo = open(path_recinfo, "a+")
    write(fid_recinfo,1807)
    write(fid_recinfo,rn) # Total Receiver Number
    write(fid_recinfo,Medium.Tn) # Total Time Sampling Tn
    write(fid_recinfo,receivers)
    close(fid_recinfo)
    #============ write Recordings ==========#
    RecLoc_Vx = (PML_Rec[2,:]-1).*WFSize.N_BDVx[1] + PML_Rec[1,:]
    RecLoc_Vz = (PML_Rec[2,:]-1).*WFSize.N_BDVz[1] + PML_Rec[1,:]
    RecLoc_P = (PML_Rec[2,:]-1).*WFSize.N_BDTpp[1] + PML_Rec[1,:]
    RecLoc_Txz = (PML_Rec[2,:]-1).*WFSize.N_BDTxz[1] + PML_Rec[1,:]
    fid_record = open(path_rec, "a+")
    for it = 1 : Tn
      Addsource!(WF1, sn, source, it)
      Onestep2D!(WF1, WF2, FD, pml)
      OptRec = [WF1.VecBDVx[RecLoc_Vx];WF1.VecBDVz[RecLoc_Vz];(WF1.VecBDTxx[RecLoc_P]+WF1.VecBDTzz[RecLoc_P])/2;WF1.VecBDTxz[RecLoc_Txz]]
      write(fid_record,OptRec)
    end
    close(fid_record)

  elseif Opt["V"] == true && Opt["P"] == true && Opt["Txz"] == false

    #========= write Recordings info =========#
    fid_recinfo = open(path_recinfo, "a+")
    write(fid_recinfo,1804)
    write(fid_recinfo,rn) # Total Receiver Number
    write(fid_recinfo,Medium.Tn) # Total Time Sampling Tn
    write(fid_recinfo,receivers)
    close(fid_recinfo)
    #============ write Recordings =============#
    RecLoc_Vx = (PML_Rec[2,:]-1).*WFSize.N_BDVx[1] + PML_Rec[1,:]
    RecLoc_Vz = (PML_Rec[2,:]-1).*WFSize.N_BDVz[1] + PML_Rec[1,:]
    RecLoc_P = (PML_Rec[2,:]-1).*WFSize.N_BDTpp[1] + PML_Rec[1,:]
    fid_record = open(path_rec, "a+")
    for it = 1 : Tn
      Addsource!(WF1, sn, source, it)
      Onestep2D!(WF1, WF2, FD, pml)
      OptRec = [WF1.VecBDVx[RecLoc_Vx];WF1.VecBDVz[RecLoc_Vz];(WF1.VecBDTxx[RecLoc_P]+WF1.VecBDTzz[RecLoc_P])/2]
      write(fid_record,OptRec)
    end
    close(fid_record)

  elseif Opt["V"] == true && Opt["P"] == false && Opt["Txz"] == true

    #========= write Recordings info =========#
    fid_recinfo = open(path_recinfo, "a+")
    write(fid_recinfo,1805)
    write(fid_recinfo,rn) # Total Receiver Number
    write(fid_recinfo,Medium.Tn) # Total Time Sampling Tn
    write(fid_recinfo,receivers)
    close(fid_recinfo)
    #============ write Recordings info ==========#
    RecLoc_Vx = (PML_Rec[2,:]-1).*WFSize.N_BDVx[1] + PML_Rec[1,:]
    RecLoc_Vz = (PML_Rec[2,:]-1).*WFSize.N_BDVz[1] + PML_Rec[1,:]
    RecLoc_Txz = (PML_Rec[2,:]-1).*WFSize.N_BDTxz[1] + PML_Rec[1,:]
    fid_record = open(path_rec, "a+")
    for it = 1 : Tn
      Addsource!(WF1, sn, source, it)
      Onestep2D!(WF1, WF2, FD, pml)
      OptRec = [WF1.VecBDVx[RecLoc_Vx];WF1.VecBDVz[RecLoc_Vz];WF1.VecBDTxz[RecLoc_Txz]]
      write(fid_record,OptRec)
    end
    close(fid_record)

  elseif Opt["V"] == false && Opt["P"] == true && Opt["Txz"] == true

    #========= write Recordings info =========#
    fid_recinfo = open(path_recinfo, "a+")
    write(fid_recinfo,1806)
    write(fid_recinfo,rn) # Total Receiver Number
    write(fid_recinfo,Medium.Tn) # Total Time Sampling Tn
    write(fid_recinfo,receivers)
    close(fid_recinfo)
    #========== write Recordings ==========#
    RecLoc_P = (PML_Rec[2,:]-1).*WFSize.N_BDTpp[1] + PML_Rec[1,:]
    RecLoc_Txz = (PML_Rec[2,:]-1).*WFSize.N_BDTxz[1] + PML_Rec[1,:]
    fid_record = open(path_rec, "a+")
    for it = 1 : Tn
      Addsource!(WF1, sn, source, it)
      Onestep2D!(WF1, WF2, FD, pml)
      OptRec = [(WF1.VecBDTxx[RecLoc_P]+WF1.VecBDTzz[RecLoc_P])/2;WF1.VecBDTxz[RecLoc_Txz]]
      write(fid_record,OptRec)
    end
    close(fid_record)

  else error("Please choose output component from V, P and Txz")
  end
end
#================================================================================================#

#======================= Simulate Single Source Seismic Recordings ==============================#
function Simulation{T,T1<:Real}(PF::T, WF1::WF2D,WF2::WF2D,WFSize::WF2DSize, Medium::Model, source::Source,
    FD::NSFDMtx, pml::dampCoef, Opt, receivers::Array{T1},rn::Int64; ext = 10, iflag = 2,
    infotxtpath = nothing, infobinpath = nothing, path_rec = "rec.bin", path_recinfo = "recinfo.bin")

  receivers = Array{Float64}(receivers)
  Rec, PML_Rec = InitRec(rn,receivers,Medium.dz,Medium.dx, ext, iflag)
  Tn = Medium.Tn
  Writeinfo(WFSize,source,PF,iflag,Medium,PML_Rec,infotxtpath,infobinpath,Opt,ext)

  if path_rec[end-2:end] != "bin" || path_recinfo[end-2:end] != "bin"
    error("Output recording data has to be in 'bin' format")
  end

  if Opt["V"] == true && Opt["P"] == false && Opt["Txz"] == false

    #========= write Recordings info =========#
    fid_recinfo = open(path_recinfo, "a+")
    write(fid_recinfo,1801)
    write(fid_recinfo,rn) # Total Receiver Number
    write(fid_recinfo,Medium.Tn) # Total Time Sampling Tn
    write(fid_recinfo,receivers)
    close(fid_recinfo)
    #========== write Recordings =========#
    RecLoc_Vx = (PML_Rec[2,:]-1).*WFSize.N_BDVx[1] + PML_Rec[1,:]
    RecLoc_Vz = (PML_Rec[2,:]-1).*WFSize.N_BDVz[1] + PML_Rec[1,:]
    fid_record = open(path_rec, "a+")
    for it = 1 : Tn
      Addsource!(WF1, source, it)
      Onestep2D!(WF1, WF2, FD, pml)
      OptRec = [WF1.VecBDVx[RecLoc_Vx];WF1.VecBDVz[RecLoc_Vz]]
      write(fid_record,OptRec)
    end
    close(fid_record)

  elseif Opt["V"] == false && Opt["P"] == true && Opt["Txz"] == false

    #========= write Recordings info =========#
    fid_recinfo = open(path_recinfo, "a+")
    write(fid_recinfo,1802)
    write(fid_recinfo,rn) # Total Receiver Number
    write(fid_recinfo,Medium.Tn) # Total Time Sampling Tn
    write(fid_recinfo,receivers)
    close(fid_recinfo)
    #=========== write Recordings ===========#
    RecLoc_P = (PML_Rec[2,:]-1).*WFSize.N_BDTpp[1] + PML_Rec[1,:]
    fid_record = open(path_rec, "a+")
    for it = 1 : Tn
      Addsource!(WF1, source, it)
      Onestep2D!(WF1, WF2, FD, pml)
      OptRec = (WF1.VecBDTxx[RecLoc_P]+WF1.VecBDTzz[RecLoc_P])/2
      write(fid_record,OptRec)
    end
    close(fid_record)

  elseif Opt["V"] == false && Opt["P"] == false && Opt["Txz"] == true

    #========= write Recordings info =========#
    fid_recinfo = open(path_recinfo, "a+")
    write(fid_recinfo,1803)
    write(fid_recinfo,rn) # Total Receiver Number
    write(fid_recinfo,Medium.Tn) # Total Time Sampling Tn
    write(fid_recinfo,receivers)
    close(fid_recinfo)
    #=========== write Recordings ===========#
    RecLoc_Txz = (PML_Rec[2,:]-1).*WFSize.N_BDTxz[1] + PML_Rec[1,:]
    fid_record = open(path_rec, "a+")
    for it = 1 : Tn
      Addsource!(WF1, source, it)
      Onestep2D!(WF1, WF2, FD, pml)
      OptRec = WF1.VecBDTxz[RecLoc_Txz]
      write(fid_record,OptRec)
    end
    close(fid_record)

  elseif Opt["V"] == true && Opt["P"] == true && Opt["Txz"] == true

    #========= write Recordings info =========#
    fid_recinfo = open(path_recinfo, "a+")
    write(fid_recinfo,1807)
    write(fid_recinfo,rn) # Total Receiver Number
    write(fid_recinfo,Medium.Tn) # Total Time Sampling Tn
    write(fid_recinfo,receivers)
    close(fid_recinfo)
    #============ write Recordings ==========#
    RecLoc_Vx = (PML_Rec[2,:]-1).*WFSize.N_BDVx[1] + PML_Rec[1,:]
    RecLoc_Vz = (PML_Rec[2,:]-1).*WFSize.N_BDVz[1] + PML_Rec[1,:]
    RecLoc_P = (PML_Rec[2,:]-1).*WFSize.N_BDTpp[1] + PML_Rec[1,:]
    RecLoc_Txz = (PML_Rec[2,:]-1).*WFSize.N_BDTxz[1] + PML_Rec[1,:]
    fid_record = open(path_rec, "a+")
    for it = 1 : Tn
      Addsource!(WF1, source, it)
      Onestep2D!(WF1, WF2, FD, pml)
      OptRec = [WF1.VecBDVx[RecLoc_Vx];WF1.VecBDVz[RecLoc_Vz];(WF1.VecBDTxx[RecLoc_P]+WF1.VecBDTzz[RecLoc_P])/2;WF1.VecBDTxz[RecLoc_Txz]]
      write(fid_record,OptRec)
    end
    close(fid_record)

  elseif Opt["V"] == true && Opt["P"] == true && Opt["Txz"] == false

    #========= write Recordings info =========#
    fid_recinfo = open(path_recinfo, "a+")
    write(fid_recinfo,1804)
    write(fid_recinfo,rn) # Total Receiver Number
    write(fid_recinfo,Medium.Tn) # Total Time Sampling Tn
    write(fid_recinfo,receivers)
    close(fid_recinfo)
    #============ write Recordings =============#
    RecLoc_Vx = (PML_Rec[2,:]-1).*WFSize.N_BDVx[1] + PML_Rec[1,:]
    RecLoc_Vz = (PML_Rec[2,:]-1).*WFSize.N_BDVz[1] + PML_Rec[1,:]
    RecLoc_P = (PML_Rec[2,:]-1).*WFSize.N_BDTpp[1] + PML_Rec[1,:]
    fid_record = open(path_rec, "a+")
    for it = 1 : Tn
      Addsource!(WF1, source, it)
      Onestep2D!(WF1, WF2, FD, pml)
      OptRec = [WF1.VecBDVx[RecLoc_Vx];WF1.VecBDVz[RecLoc_Vz];(WF1.VecBDTxx[RecLoc_P]+WF1.VecBDTzz[RecLoc_P])/2]
      write(fid_record,OptRec)
    end
    close(fid_record)

  elseif Opt["V"] == true && Opt["P"] == false && Opt["Txz"] == true

    #========= write Recordings info =========#
    fid_recinfo = open(path_recinfo, "a+")
    write(fid_recinfo,1805)
    write(fid_recinfo,rn) # Total Receiver Number
    write(fid_recinfo,Medium.Tn) # Total Time Sampling Tn
    write(fid_recinfo,receivers)
    close(fid_recinfo)
    #============ write Recordings info ==========#
    RecLoc_Vx = (PML_Rec[2,:]-1).*WFSize.N_BDVx[1] + PML_Rec[1,:]
    RecLoc_Vz = (PML_Rec[2,:]-1).*WFSize.N_BDVz[1] + PML_Rec[1,:]
    RecLoc_Txz = (PML_Rec[2,:]-1).*WFSize.N_BDTxz[1] + PML_Rec[1,:]
    fid_record = open(path_rec, "a+")
    for it = 1 : Tn
      Addsource!(WF1, source, it)
      Onestep2D!(WF1, WF2, FD, pml)
      OptRec = [WF1.VecBDVx[RecLoc_Vx];WF1.VecBDVz[RecLoc_Vz];WF1.VecBDTxz[RecLoc_Txz]]
      write(fid_record,OptRec)
    end
    close(fid_record)

  elseif Opt["V"] == false && Opt["P"] == true && Opt["Txz"] == true

    #========= write Recordings info =========#
    fid_recinfo = open(path_recinfo, "a+")
    write(fid_recinfo,1806)
    write(fid_recinfo,rn) # Total Receiver Number
    write(fid_recinfo,Medium.Tn) # Total Time Sampling Tn
    write(fid_recinfo,receivers)
    close(fid_recinfo)
    #========== write Recordings ==========#
    RecLoc_P = (PML_Rec[2,:]-1).*WFSize.N_BDTpp[1] + PML_Rec[1,:]
    RecLoc_Txz = (PML_Rec[2,:]-1).*WFSize.N_BDTxz[1] + PML_Rec[1,:]
    fid_record = open(path_rec, "a+")
    for it = 1 : Tn
      Addsource!(WF1, source, it)
      Onestep2D!(WF1, WF2, FD, pml)
      OptRec = [(WF1.VecBDTxx[RecLoc_P]+WF1.VecBDTzz[RecLoc_P])/2;WF1.VecBDTxz[RecLoc_Txz]]
      write(fid_record,OptRec)
    end
    close(fid_record)

  else error("Please choose output component from V, P and Txz")
  end
end
#================================================================================================#

#======================= Simulate Multiple Source Seismic Wavefields ==============================#
function Simulation{T<:Real}(PF::T, WF1::WF2D,WF2::WF2D,WFSize::WF2DSize, Medium::Model, source::Array{Source},
  FD::NSFDMtx, pml::dampCoef, sn::Int64, Opt, Slices::Union{UnitRange{Int64},StepRange{Int64,Int64},Array{Int64,1},Int64};
  ext = 10, iflag = 2, infotxtpath = nothing, infobinpath = nothing, path_wf = "wf.bin", path_wfinfo = "wfinfo.bin")

Tn = Medium.Tn
Writeinfo(WFSize,source,sn,PF,iflag,Medium,Slices,infotxtpath,infobinpath,Opt,ext)

if path_wf[end-2:end] != "bin" || path_wfinfo[end-2:end] != "bin"
  error("Output wavefield data has to be in 'bin' format")
end

  if Opt["V"] == true && Opt["P"] == false && Opt["Txz"] == false

    #==== write Wavefields info ===========#
    fid_wfinfo = open(path_wfinfo, "a+")
    write(fid_wfinfo,1801)
    write(fid_wfinfo,WFSize.N_BDTpp)
    write(fid_wfinfo,length(Slices))
    write(fid_wfinfo,Slices.*Medium.dt)
    close(fid_wfinfo)
    #============  write Wavefields ===========#
    fid_wf = open(path_wf, "a+")
    for it = 1 : Tn
      Addsource!(WF1, sn, source, it)
      Onestep2D!(WF1, WF2, FD, pml)
      if it in Slices
        OptWF = [WF1.VecBDVx;WF1.VecBDVz]
        write(fid_wf,OptWF)
      end
    end
    close(fid_wf)

  elseif Opt["V"] == false && Opt["P"] == true && Opt["Txz"] == false

    #========= write Wavefields info ===========#
    fid_wfinfo = open(path_wfinfo, "a+")
    write(fid_wfinfo,1802)
    write(fid_wfinfo,WFSize.N_BDTpp)
    write(fid_wfinfo,length(Slices))
    write(fid_wfinfo,Slices.*Medium.dt)
    close(fid_wfinfo)
    #============= write Wavefields ============#
    fid_wf = open(path_wf, "a+")
    for it = 1 : Tn
      Addsource!(WF1, sn, source, it)
      Onestep2D!(WF1, WF2, FD, pml)
      if it in Slices
        OptWF = (WF1.VecBDTxx+WF1.VecBDTzz)/2
        write(fid_wf,OptWF)
      end
    end
    close(fid_wf)

  elseif Opt["V"] == false && Opt["P"] == false && Opt["Txz"] == true

    #========= write Wavefields info ===========#
    fid_wfinfo = open(path_wfinfo, "a+")
    write(fid_wfinfo,1803)
    write(fid_wfinfo,WFSize.N_BDTpp)
    write(fid_wfinfo,length(Slices))
    write(fid_wfinfo,Slices.*Medium.dt)
    close(fid_wfinfo)
    #==========  write Wavefields ===========#
    fid_wf = open(path_wf, "a+")
    for it = 1 : Tn
      Addsource!(WF1, sn, source, it)
      Onestep2D!(WF1, WF2, FD, pml)
      if it in Slices
        OptWF = WF1.VecBDTxz
        write(fid_wf,OptWF)
      end
    end
    close(fid_wf)

  elseif Opt["V"] == true && Opt["P"] == true && Opt["Txz"] == true

    #========= write Wavefields info ===========#
    fid_wfinfo = open(path_wfinfo, "a+")
    write(fid_wfinfo,1807)
    write(fid_wfinfo,WFSize.N_BDTpp)
    write(fid_wfinfo,length(Slices))
    write(fid_wfinfo,Slices.*Medium.dt)
    close(fid_wfinfo)
    #============ write Wavefields ============#
    fid_wf = open(path_wf, "a+")
    for it = 1 : Tn
      Addsource!(WF1, sn, source, it)
      Onestep2D!(WF1, WF2, FD, pml)
      if it in Slices
        OptWF = [WF1.VecBDVx;WF1.VecBDVz;(WF1.VecBDTxx+WF1.VecBDTzz)/2;WF1.VecBDTxz]
        write(fid_wf,OptWF)
      end
    end
    close(fid_wf)

  elseif Opt["V"] == true && Opt["P"] == true && Opt["Txz"] == false

    #========= write Wavefields info===========#
    fid_wfinfo = open(path_wfinfo, "a+")
    write(fid_wfinfo,1804)
    write(fid_wfinfo,WFSize.N_BDTpp)
    write(fid_wfinfo,length(Slices))
    write(fid_wfinfo,Slices.*Medium.dt)
    close(fid_wfinfo)
    #============  write Wavefields ============#
    fid_wf = open(path_wf, "a+")
    for it = 1 : Tn
      Addsource!(WF1, sn, source, it)
      Onestep2D!(WF1, WF2, FD, pml)
      if it in Slices
        OptWF = [WF1.VecBDVx;WF1.VecBDVz;(WF1.VecBDTxx+WF1.VecBDTzz)/2]
        write(fid_wf,OptWF)
      end
    end
    close(fid_wf)

  elseif Opt["V"] == true && Opt["P"] == false && Opt["Txz"] == true

    #========= write Wavefields info ===========#
    fid_wfinfo = open(path_wfinfo, "a+")
    write(fid_wfinfo,1805)
    write(fid_wfinfo,WFSize.N_BDTpp)
    write(fid_wfinfo,length(Slices))
    write(fid_wfinfo,Slices.*Medium.dt)
    close(fid_wfinfo)
    #===========  write Wavefields =============#
    fid_wf = open(path_wf, "a+")
    for it = 1 : Tn
      Addsource!(WF1, sn, source, it)
      Onestep2D!(WF1, WF2, FD, pml)
      if it in Slices
        OptWF = [WF1.VecBDVx;WF1.VecBDVz;WF1.VecBDTxz]
        write(fid_wf,OptWF)
      end
    end
    close(fid_wf)

  elseif Opt["V"] == false && Opt["P"] == true && Opt["Txz"] == true

    #========= write Wavefields info ===========#
    fid_wfinfo = open(path_wfinfo, "a+")
    write(fid_wfinfo,1806)
    write(fid_wfinfo,WFSize.N_BDTpp)
    write(fid_wfinfo,length(Slices))
    write(fid_wfinfo,Slices.*Medium.dt)
    close(fid_wfinfo)
    #===========  write Wavefields ============#
    fid_wf = open(path_wf, "a+")
    for it = 1 : Tn
      Addsource!(WF1, sn, source, it)
      Onestep2D!(WF1, WF2, FD, pml)
      if it in Slices
        OptWF = [(WF1.VecBDTxx+WF1.VecBDTzz)/2;WF1.VecBDTxz]
        write(fid_wf,OptWF)
      end
    end
    close(fid_wf)

  else error("Please choose output component from V, P and Txz")
  end
end
#================================================================================================#

#======================= Simulate Single Source Seismic Wavefields ==============================#
function Simulation{T<:Real}(PF::T, WF1::WF2D,WF2::WF2D,WFSize::WF2DSize, Medium::Model, source::Source,
  FD::NSFDMtx, pml::dampCoef, Opt,
  Slices::Union{UnitRange{Int64},StepRange{Int64,Int64},Array{Int64,1}, Int64}; ext = 10, iflag = 2,
  infotxtpath = nothing, infobinpath = nothing, path_wf = "wf.bin", path_wfinfo = "wfinfo.bin")

  Tn = Medium.Tn
  Writeinfo(WFSize,source,PF,iflag,Medium,Slices,infotxtpath,infobinpath,Opt,ext)

  if path_wf[end-2:end] != "bin" || path_wfinfo[end-2:end] != "bin"
  error("Output wavefield data has to be in 'bin' format")
  end

    if Opt["V"] == true && Opt["P"] == false && Opt["Txz"] == false

      #========= write Wavefields info ===========#
      fid_wfinfo = open(path_wfinfo, "a+")
      write(fid_wfinfo,1801)
      write(fid_wfinfo,WFSize.N_BDTpp)
      write(fid_wfinfo,length(Slices))
      write(fid_wfinfo,Slices.*Medium.dt)
      close(fid_wfinfo)
      #============  write Wavefields ===========#
      fid_wf = open(path_wf, "a+")
      for it = 1 : Tn
        Addsource!(WF1, source, it)
        Onestep2D!(WF1, WF2, FD, pml)
        if it in Slices
          OptWF = [WF1.VecBDVx;WF1.VecBDVz]
          write(fid_wf,OptWF)
        end
      end
      close(fid_wf)

    elseif Opt["V"] == false && Opt["P"] == true && Opt["Txz"] == false

      #========= write Wavefields info ===========#
      fid_wfinfo = open(path_wfinfo, "a+")
      write(fid_wfinfo,1802)
      write(fid_wfinfo,WFSize.N_BDTpp)
      write(fid_wfinfo,length(Slices))
      write(fid_wfinfo,Slices.*Medium.dt)
      close(fid_wfinfo)
      #============= write Wavefields ============#
      fid_wf = open(path_wf, "a+")
      for it = 1 : Tn
        Addsource!(WF1, source, it)
        Onestep2D!(WF1, WF2, FD, pml)
        if it in Slices
          OptWF = (WF1.VecBDTxx+WF1.VecBDTzz)/2
          write(fid_wf,OptWF)
        end
      end
      close(fid_wf)

    elseif Opt["V"] == false && Opt["P"] == false && Opt["Txz"] == true

      #========= write Wavefields info ===========#
      fid_wfinfo = open(path_wfinfo, "a+")
      write(fid_wfinfo,1803)
      write(fid_wfinfo,WFSize.N_BDTpp)
      write(fid_wfinfo,length(Slices))
      write(fid_wfinfo,Slices.*Medium.dt)
      close(fid_wfinfo)
      #==========  write Wavefields ===========#
      fid_wf = open(path_wf, "a+")
      for it = 1 : Tn
        Addsource!(WF1, source, it)
        Onestep2D!(WF1, WF2, FD, pml)
        if it in Slices
          OptWF = WF1.VecBDTxz
          write(fid_wf,OptWF)
        end
      end
      close(fid_wf)

    elseif Opt["V"] == true && Opt["P"] == true && Opt["Txz"] == true

      #========= write Wavefields info ===========#
      fid_wfinfo = open(path_wfinfo, "a+")
      write(fid_wfinfo,1807)
      write(fid_wfinfo,WFSize.N_BDTpp)
      write(fid_wfinfo,length(Slices))
      write(fid_wfinfo,Slices.*Medium.dt)
      close(fid_wfinfo)
      #============ write Wavefields ============#
      fid_wf = open(path_wf, "a+")
      for it = 1 : Tn
        Addsource!(WF1, source, it)
        Onestep2D!(WF1, WF2, FD, pml)
        if it in Slices
          OptWF = [WF1.VecBDVx;WF1.VecBDVz;(WF1.VecBDTxx+WF1.VecBDTzz)/2;WF1.VecBDTxz]
          write(fid_wf,OptWF)
        end
      end
      close(fid_wf)

    elseif Opt["V"] == true && Opt["P"] == true && Opt["Txz"] == false

      #========= write Wavefields info===========#
      fid_wfinfo = open(path_wfinfo, "a+")
      write(fid_wfinfo,1804)
      write(fid_wfinfo,WFSize.N_BDTpp)
      write(fid_wfinfo,length(Slices))
      write(fid_wfinfo,Slices.*Medium.dt)
      close(fid_wfinfo)
      #============  write Wavefields ============#
      fid_wf = open(path_wf, "a+")
      for it = 1 : Tn
        Addsource!(WF1, source, it)
        Onestep2D!(WF1, WF2, FD, pml)
        if it in Slices
          OptWF = [WF1.VecBDVx;WF1.VecBDVz;(WF1.VecBDTxx+WF1.VecBDTzz)/2]
          write(fid_wf,OptWF)
        end
      end
      close(fid_wf)

    elseif Opt["V"] == true && Opt["P"] == false && Opt["Txz"] == true

      #========= write Wavefields info ===========#
      fid_wfinfo = open(path_wfinfo, "a+")
      write(fid_wfinfo,1805)
      write(fid_wfinfo,WFSize.N_BDTpp)
      write(fid_wfinfo,length(Slices))
      write(fid_wfinfo,Slices.*Medium.dt)
      close(fid_wfinfo)
      #===========  write Wavefields =============#
      fid_wf = open(path_wf, "a+")
      for it = 1 : Tn
        Addsource!(WF1, source, it)
        Onestep2D!(WF1, WF2, FD, pml)
        if it in Slices
          OptWF = [WF1.VecBDVx;WF1.VecBDVz;WF1.VecBDTxz]
          write(fid_wf,OptWF)
        end
      end
      close(fid_wf)

    elseif Opt["V"] == false && Opt["P"] == true && Opt["Txz"] == true

      #========= write Wavefields info ===========#
      fid_wfinfo = open(path_wfinfo, "a+")
      write(fid_wfinfo,1806)
      write(fid_wfinfo,WFSize.N_BDTpp)
      write(fid_wfinfo,length(Slices))
      write(fid_wfinfo,Slices.*Medium.dt)
      close(fid_wfinfo)
      #===========  write Wavefields ============#
      fid_wf = open(path_wf, "a+")
      for it = 1 : Tn
        Addsource!(WF1, source, it)
        Onestep2D!(WF1, WF2, FD, pml)
        if it in Slices
          OptWF = [(WF1.VecBDTxx+WF1.VecBDTzz)/2;WF1.VecBDTxz]
          write(fid_wf,OptWF)
        end
      end
      close(fid_wf)

    else error("Please choose output component from V, P and Txz")
    end
  end
#================================================================================================#

#======================= Simulate Multiple Source Seismic Wavefields for Model Order Reduction ==============================#
function Simulation{T<:Real}(PF::T, WF1::WF2D,WF2::WF2D,WFSize::WF2DSize, Medium::Model, source::Array{Source},
  FD::SFDMtx, sn::Int64, Slices::Union{UnitRange{Int64},StepRange{Int64,Int64},Array{Int64,1},Int64};
  ext = 10, iflag = 2, infotxtpath = nothing, infobinpath = nothing, path_wf = "wf.bin", path_wfinfo = nothing)

    Tn = Medium.Tn
    Writeinfo(WFSize,source,sn,PF,iflag,Medium,Slices,infotxtpath,infobinpath,ext)
    if path_wfinfo != nothing
      if path_wf[end-2:end] != "bin" || path_wfinfo[end-2:end] != "bin"
        error("Output wavefield data has to be in 'bin' format")
      end
    #========= write Wavefields info ===========#
    fid_wfinfo = open(path_wfinfo, "a+")
    write(fid_wfinfo,1808)
    write(fid_wfinfo,WFSize.N_BDTpp)
    write(fid_wfinfo,length(Slices))
    write(fid_wfinfo,Slices.*Medium.dt)
    close(fid_wfinfo)
  end
    #============ write Wavefields ============#
    fid_wf = open(path_wf, "a+")
    for it = 1 : Tn
      AddSLsource!(WF1, sn, source, it)
      Onestep2D!(WF1, WF2, FD)
      if it in Slices
        Vx = [WF1.VecBDVxx ; WF1.VecBDVxz]
        Vz = [WF1.VecBDVzx ; WF1.VecBDVzz]
        Txx = [WF1.VecBDTxxx ; WF1.VecBDTxxz]
        Tzz = [WF1.VecBDTzzx ; WF1.VecBDTzzz]
        Txz = [WF1.VecBDTxzx ; WF1.VecBDTxzz]
        OptWF = [Vx;Vz;Txx;Tzz;Txz]
        write(fid_wf,OptWF)
      end
    end
    close(fid_wf)

end
#================================================================================================#

#======================= Simulate Single Source Seismic Wavefields for Model Order Reduction ==============================#
function Simulation{T<:Real}(PF::T, WF1::WF2D, WF2::WF2D, WFSize::WF2DSize, Medium::Model, source::Source,
  FD::SFDMtx, Slices::Union{UnitRange{Int64},StepRange{Int64,Int64},Array{Int64,1}, Int64}; ext = 10, iflag = 2,
  infotxtpath = nothing, infobinpath = nothing, path_wf = "wf.bin", path_wfinfo = nothing)

  Tn = Medium.Tn
  Writeinfo(WFSize,source,PF,iflag,Medium,Slices,infotxtpath,infobinpath,ext)
      if path_wfinfo != nothing
        if path_wf[end-2:end] != "bin" || path_wfinfo[end-2:end] != "bin"
          error("Output wavefield data has to be in 'bin' format")
        end
      #========= write Wavefields info ===========#
        fid_wfinfo = open(path_wfinfo, "a+")
        write(fid_wfinfo,1808)
        write(fid_wfinfo,WFSize.N_BDTpp)
        write(fid_wfinfo,length(Slices))
        write(fid_wfinfo,Slices.*Medium.dt)
        close(fid_wfinfo)
      end
      #============ write Wavefields ============#
      fid_wf = open(path_wf, "a+")
      for it = 1 : Tn
        AddSLsource!(WF1, source, it)
        Onestep2D!(WF1, WF2, FD)
        if it in Slices
          Vx = [WF1.VecBDVxx ; WF1.VecBDVxz]
          Vz = [WF1.VecBDVzx ; WF1.VecBDVzz]
          Txx = [WF1.VecBDTxxx ; WF1.VecBDTxxz]
          Tzz = [WF1.VecBDTzzx ; WF1.VecBDTzzz]
          Txz = [WF1.VecBDTxzx ; WF1.VecBDTxzz]
          OptWF = [Vx;Vz;Txx;Tzz;Txz]
          write(fid_wf,OptWF)
        end
      end
      close(fid_wf)

  end
#================================================================================================#

#================= Simulate Multiple Source Seismic Recordings and Wavefields ===================#
function Simulation{T,T1<:Real}(PF::T, WF1::WF2D,WF2::WF2D,WFSize::WF2DSize, Medium::Model, source::Array{Source},
   FD::NSFDMtx, pml::dampCoef, sn::Int64, Opt, receivers::Array{T1},rn::Int64,
   Slices::Union{UnitRange{Int64},StepRange{Int64,Int64},Array{Int64,1}, Int64}; ext = 10, iflag = 2,
   infotxtpath = nothing, infobinpath = nothing, path_wf = "wf.bin", path_wfinfo = "wfinfo.bin",
   path_rec = "rec.bin", path_recinfo = "recinfo.bin")

   receivers = Array{Float64}(receivers)
   Rec, PML_Rec = InitRec(rn,receivers,Medium.dz,Medium.dx, ext, iflag)
  Tn = Medium.Tn
  Writeinfo(WFSize,source,sn,PF,iflag,Medium,PML_Rec,Slices,infotxtpath,infobinpath, Opt,ext)

  if path_rec[end-2:end] != "bin" || path_recinfo[end-2:end] != "bin"
    error("Output recording data has to be in 'bin' format")
  end

  if path_wf[end-2:end] != "bin" || path_wfinfo[end-2:end] != "bin"
    error("Output wavefield data has to be in 'bin' format")
  end

  if Opt["V"] == true && Opt["P"] == false && Opt["Txz"] == false

    #========= write Recordings.bin =========#
    fid_recinfo = open(path_recinfo, "a+")
    write(fid_recinfo,1801)
    write(fid_recinfo,rn) # Total Receiver Number
    write(fid_recinfo,Medium.Tn) # Total Time Sampling Tn
    write(fid_recinfo,receivers)
    close(fid_recinfo)
    #==========================================#

    #========= write Wavefields info ===========#
    fid_wfinfo = open(path_wfinfo, "a+")
    write(fid_wfinfo,1801)
    write(fid_wfinfo,WFSize.N_BDTpp)
    write(fid_wfinfo,length(Slices))
    write(fid_wfinfo,Slices.*Medium.dt)
    close(fid_wfinfo)
    #==========================================#
    RecLoc_Vx = (PML_Rec[2,:]-1).*WFSize.N_BDVx[1] + PML_Rec[1,:]
    RecLoc_Vz = (PML_Rec[2,:]-1).*WFSize.N_BDVz[1] + PML_Rec[1,:]
    fid_record = open(path_rec, "a+")
    fid_wf = open(path_wf, "a+")
    for it = 1 : Tn
      Addsource!(WF1, sn, source, it)
      Onestep2D!(WF1, WF2, FD, pml)
      OptRec = [WF1.VecBDVx[RecLoc_Vx];WF1.VecBDVz[RecLoc_Vz]]
      write(fid_record,OptRec)
      if it in Slices
        OptWF = [WF1.VecBDVx;WF1.VecBDVz]
        write(fid_wf,OptWF)
      end
    end
    close(fid_record)
    close(fid_wf)

  elseif Opt["V"] == false && Opt["P"] == true && Opt["Txz"] == false

    #========= write Recordings info =========#
    fid_recinfo = open(path_recinfo, "a+")
    write(fid_recinfo,1802)
    write(fid_recinfo,rn) # Total Receiver Number
    write(fid_recinfo,Medium.Tn) # Total Time Sampling Tn
    write(fid_recinfo,receivers)
    close(fid_recinfo)
    #==========================================#

    #========= write Wavefields info ===========#
    fid_wfinfo = open(path_wfinfo, "a+")
    write(fid_wfinfo,1802)
    write(fid_wfinfo,WFSize.N_BDTpp)
    write(fid_wfinfo,length(Slices))
    write(fid_wfinfo,Slices.*Medium.dt)
    close(fid_wfinfo)
    #==========================================#
    RecLoc_P = (PML_Rec[2,:]-1).*WFSize.N_BDTpp[1] + PML_Rec[1,:]
    fid_record = open(path_rec,"a+")
    fid_wf = open(path_wf,"a+")
    for it = 1 : Tn
      Addsource!(WF1, sn, source, it)
      Onestep2D!(WF1, WF2, FD, pml)
      OptRec = (WF1.VecBDTxx[RecLoc_P]+WF1.VecBDTzz[RecLoc_P])/2
      write(fid_record,OptRec)
      if it in Slices
        OptWF = (WF1.VecBDTxx+WF1.VecBDTzz)/2
        write(fid_wf,OptWF)
      end
    end
    close(fid_record)
    close(fid_wf)

  elseif Opt["V"] == false && Opt["P"] == false && Opt["Txz"] == true

    #========= write Recordings info =========#
    fid_recinfo = open(path_recinfo, "a+")
    write(fid_recinfo,1803)
    write(fid_recinfo,rn) # Total Receiver Number
    write(fid_recinfo,Medium.Tn) # Total Time Sampling Tn
    write(fid_recinfo,receivers)
    close(fid_recinfo)
    #==========================================#

    #========= write Wavefields info ===========#
    fid_wfinfo = open(path_wfinfo, "a+")
    write(fid_wfinfo,1803)
    write(fid_wfinfo,WFSize.N_BDTpp)
    write(fid_wfinfo,length(Slices))
    write(fid_wfinfo,Slices.*Medium.dt)
    close(fid_wfinfo)
    #==========================================#
    RecLoc_Txz = (PML_Rec[2,:]-1).*WFSize.N_BDTxz[1] + PML_Rec[1,:]
    fid_record = open(path_rec,"a+")
    fid_wf = open(path_wf,"a+")
    for it = 1 : Tn
      Addsource!(WF1, sn, source, it)
      Onestep2D!(WF1, WF2, FD, pml)
      OptRec = WF1.VecBDTxz[RecLoc_Txz]
      write(fid_record,OptRec)
      if it in Slices
        OptWF = WF1.VecBDTxz
        write(fid_wf,OptWF)
      end
    end
    close(fid_record)
    close(fid_wf)

  elseif Opt["V"] == true && Opt["P"] == true && Opt["Txz"] == true

    #========= write Recordings info =========#
    fid_recinfo = open(path_recinfo, "a+")
    write(fid_recinfo,1807)
    write(fid_recinfo,rn) # Total Receiver Number
    write(fid_recinfo,Medium.Tn) # Total Time Sampling Tn
    write(fid_recinfo,receivers)
    close(fid_recinfo)
    #==========================================#

    #========= write Wavefields info ===========#
    fid_wfinfo = open(path_wfinfo, "a+")
    write(fid_wfinfo,1807)
    write(fid_wfinfo,WFSize.N_BDTpp)
    write(fid_wfinfo,length(Slices))
    write(fid_wfinfo,Slices.*Medium.dt)
    close(fid_wfinfo)
    #==========================================#
    RecLoc_Vx = (PML_Rec[2,:]-1).*WFSize.N_BDVx[1] + PML_Rec[1,:]
    RecLoc_Vz = (PML_Rec[2,:]-1).*WFSize.N_BDVz[1] + PML_Rec[1,:]
    RecLoc_P = (PML_Rec[2,:]-1).*WFSize.N_BDTpp[1] + PML_Rec[1,:]
    RecLoc_Txz = (PML_Rec[2,:]-1).*WFSize.N_BDTxz[1] + PML_Rec[1,:]
    fid_record = open(path_rec,"a+")
    fid_wf = open(path_wf,"a+")
    for it = 1 : Tn
      Addsource!(WF1, sn, source, it)
      Onestep2D!(WF1, WF2, FD, pml)
      OptRec = [WF1.VecBDVx[RecLoc_Vx];WF1.VecBDVz[RecLoc_Vz];(WF1.VecBDTxx[RecLoc_P]+WF1.VecBDTzz[RecLoc_P])/2;WF1.VecBDTxz[RecLoc_Txz]]
      write(fid_record,OptRec)
      if it in Slices
        OptWF = [WF1.VecBDVx;WF1.VecBDVz;(WF1.VecBDTxx+WF1.VecBDTzz)/2;WF1.VecBDTxz]
        write(fid_wf,OptWF)
      end
    end
    close(fid_record)
    close(fid_wf)

  elseif Opt["V"] == true && Opt["P"] == true && Opt["Txz"] == false

    #========= write Recordings info =========#
    fid_recinfo = open(path_recinfo, "a+")
    write(fid_recinfo,1804)
    write(fid_recinfo,rn) # Total Receiver Number
    write(fid_recinfo,Medium.Tn) # Total Time Sampling Tn
    write(fid_recinfo,receivers)
    close(fid_recinfo)
    #==========================================#

    #========= write Wavefields info===========#
    fid_wfinfo = open(path_wfinfo, "a+")
    write(fid_wfinfo,1804)
    write(fid_wfinfo,WFSize.N_BDTpp)
    write(fid_wfinfo,length(Slices))
    write(fid_wfinfo,Slices.*Medium.dt)
    close(fid_wfinfo)
    #==========================================#
    RecLoc_Vx = (PML_Rec[2,:]-1).*WFSize.N_BDVx[1] + PML_Rec[1,:]
    RecLoc_Vz = (PML_Rec[2,:]-1).*WFSize.N_BDVz[1] + PML_Rec[1,:]
    RecLoc_P = (PML_Rec[2,:]-1).*WFSize.N_BDTpp[1] + PML_Rec[1,:]
    fid_record = open(path_rec,"a+")
    fid_wf = open(path_wf,"a+")
    for it = 1 : Tn
      Addsource!(WF1, sn, source, it)
      Onestep2D!(WF1, WF2, FD, pml)
      OptRec = [WF1.VecBDVx[RecLoc_Vx];WF1.VecBDVz[RecLoc_Vz];(WF1.VecBDTxx[RecLoc_P]+WF1.VecBDTzz[RecLoc_P])/2]
      write(fid_record,OptRec)
      if it in Slices
        OptWF = [WF1.VecBDVx;WF1.VecBDVz;(WF1.VecBDTxx+WF1.VecBDTzz)/2]
        write(fid_wf,OptWF)
      end
    end
    close(fid_record)
    close(fid_wf)

  elseif Opt["V"] == true && Opt["P"] == false && Opt["Txz"] == true

    #========= write Recordings info =========#
    fid_recinfo = open(path_recinfo, "a+")
    write(fid_recinfo,1805)
    write(fid_recinfo,rn) # Total Receiver Number
    write(fid_recinfo,Medium.Tn) # Total Time Sampling Tn
    write(fid_recinfo,receivers)
    close(fid_recinfo)
    #==========================================#

    #========= write Wavefields info ===========#
    fid_wfinfo = open(path_wfinfo, "a+")
    write(fid_wfinfo,1805)
    write(fid_wfinfo,WFSize.N_BDTpp)
    write(fid_wfinfo,length(Slices))
    write(fid_wfinfo,Slices.*Medium.dt)
    close(fid_wfinfo)
    #==========================================#
    RecLoc_Vx = (PML_Rec[2,:]-1).*WFSize.N_BDVx[1] + PML_Rec[1,:]
    RecLoc_Vz = (PML_Rec[2,:]-1).*WFSize.N_BDVz[1] + PML_Rec[1,:]
    RecLoc_Txz = (PML_Rec[2,:]-1).*WFSize.N_BDTxz[1] + PML_Rec[1,:]
    fid_record = open(path_rec,"a+")
    fid_wf = open(path_wf,"a+")
    for it = 1 : Tn
      Addsource!(WF1, sn, source, it)
      Onestep2D!(WF1, WF2, FD, pml)
      OptRec = [WF1.VecBDVx[RecLoc_Vx];WF1.VecBDVz[RecLoc_Vz];WF1.VecBDTxz[RecLoc_Txz]]
      write(fid_record,OptRec)
      if it in Slices
        OptWF = [WF1.VecBDVx;WF1.VecBDVz;WF1.VecBDTxz]
        write(fid_wf,OptWF)
      end
    end
    close(fid_record)
    close(fid_wf)

  elseif Opt["V"] == false && Opt["P"] == true && Opt["Txz"] == true

    #========= write Recordings info =========#
    fid_recinfo = open(path_recinfo, "a+")
    write(fid_recinfo,1806)
    write(fid_recinfo,rn) # Total Receiver Number
    write(fid_recinfo,Medium.Tn) # Total Time Sampling Tn
    write(fid_recinfo,receivers)
    close(fid_recinfo)
    #==========================================#

    #========= write Wavefields info ===========#
    fid_wfinfo = open(path_wfinfo, "a+")
    write(fid_wfinfo,1806)
    write(fid_wfinfo,WFSize.N_BDTpp)
    write(fid_wfinfo,length(Slices))
    write(fid_wfinfo,Slices.*Medium.dt)
    close(fid_wfinfo)
    #==========================================#
    RecLoc_P = (PML_Rec[2,:]-1).*WFSize.N_BDTpp[1] + PML_Rec[1,:]
    RecLoc_Txz = (PML_Rec[2,:]-1).*WFSize.N_BDTxz[1] + PML_Rec[1,:]
    fid_record = open(path_rec,"a+")
    fid_wf = open(path_wf,"a+")
    for it = 1 : Tn
      Addsource!(WF1, sn, source, it)
      Onestep2D!(WF1, WF2, FD, pml)
      OptRec = [(WF1.VecBDTxx[RecLoc_P]+WF1.VecBDTzz[RecLoc_P])/2;WF1.VecBDTxz[RecLoc_Txz]]
      write(fid_record,OptRec)
      if it in Slices
        OptWF = [(WF1.VecBDTxx+WF1.VecBDTzz)/2;WF1.VecBDTxz]
        write(fid_wf,OptWF)
      end
    end
    close(fid_record)
    close(fid_wf)

  else error("Please choose output component from V, P and Txz")
  end
end
#================================================================================================#

#================= Simulate Single Source Seismic Recordings and Wavefields ===================#
function Simulation{T,T1<:Real}(PF::T, WF1::WF2D,WF2::WF2D,WFSize::WF2DSize, Medium::Model, source::Source,
   FD::NSFDMtx, pml::dampCoef, Opt, receivers::Array{T1},rn::Int64,
   Slices::Union{UnitRange{Int64},StepRange{Int64,Int64},Array{Int64,1}, Int64}; ext = 10, iflag = 2,
   infotxtpath = nothing, infobinpath = nothing, path_wf = "wf.bin", path_wfinfo = "wfinfo.bin",
   path_rec = "rec.bin", path_recinfo = "recinfo.bin")

   receivers = Array{Float64}(receivers)
   Rec, PML_Rec = InitRec(rn,receivers,Medium.dz,Medium.dx, ext, iflag)
   Tn = Medium.Tn
   Writeinfo(WFSize,source,sn,PF,iflag,Medium,PML_Rec,Slices,infotxtpath,infobinpath, Opt,ext)

   if path_rec[end-2:end] != "bin" || path_recinfo[end-2:end] != "bin"
    error("Output recording data has to be in 'bin' format")
  end

  if path_wf[end-2:end] != "bin" || path_wfinfo[end-2:end] != "bin"
    error("Output wavefield data has to be in 'bin' format")
  end

   if Opt["V"] == true && Opt["P"] == false && Opt["Txz"] == false

     #========= write Recordings.bin =========#
     fid_recinfo = open(path_recinfo, "a+")
     write(fid_recinfo,1801)
     write(fid_recinfo,rn) # Total Receiver Number
     write(fid_recinfo,Medium.Tn) # Total Time Sampling Tn
     write(fid_recinfo,receivers)
     close(fid_recinfo)
     #==========================================#

     #========= write Wavefields info ===========#
     fid_wfinfo = open(path_wfinfo, "a+")
     write(fid_wfinfo,1801)
     write(fid_wfinfo,WFSize.N_BDTpp)
     write(fid_wfinfo,length(Slices))
     write(fid_wfinfo,Slices.*Medium.dt)
     close(fid_wfinfo)
     #==========================================#
     RecLoc_Vx = (PML_Rec[2,:]-1).*WFSize.N_BDVx[1] + PML_Rec[1,:]
     RecLoc_Vz = (PML_Rec[2,:]-1).*WFSize.N_BDVz[1] + PML_Rec[1,:]
     fid_record = open(path_rec, "a+")
     fid_wf = open(path_wf, "a+")
     for it = 1 : Tn
       Addsource!(WF1, source, it)
       Onestep2D!(WF1, WF2, FD, pml)
       OptRec = [WF1.VecBDVx[RecLoc_Vx];WF1.VecBDVz[RecLoc_Vz]]
       write(fid_record,OptRec)
       if it in Slices
         OptWF = [WF1.VecBDVx;WF1.VecBDVz]
         write(fid_wf,OptWF)
       end
     end
     close(fid_record)
     close(fid_wf)

   elseif Opt["V"] == false && Opt["P"] == true && Opt["Txz"] == false

     #========= write Recordings info =========#
     fid_recinfo = open(path_recinfo, "a+")
     write(fid_recinfo,1802)
     write(fid_recinfo,rn) # Total Receiver Number
     write(fid_recinfo,Medium.Tn) # Total Time Sampling Tn
     write(fid_recinfo,receivers)
     close(fid_recinfo)
     #==========================================#

     #========= write Wavefields info ===========#
     fid_wfinfo = open(path_wfinfo, "a+")
     write(fid_wfinfo,1802)
     write(fid_wfinfo,WFSize.N_BDTpp)
     write(fid_wfinfo,length(Slices))
     write(fid_wfinfo,Slices.*Medium.dt)
     close(fid_wfinfo)
     #==========================================#
     RecLoc_P = (PML_Rec[2,:]-1).*WFSize.N_BDTpp[1] + PML_Rec[1,:]
     fid_record = open(path_rec,"a+")
     fid_wf = open(path_wf,"a+")
     for it = 1 : Tn
       Addsource!(WF1, source, it)
       Onestep2D!(WF1, WF2, FD, pml)
       OptRec = (WF1.VecBDTxx[RecLoc_P]+WF1.VecBDTzz[RecLoc_P])/2
       write(fid_record,OptRec)
       if it in Slices
         OptWF = (WF1.VecBDTxx+WF1.VecBDTzz)/2
         write(fid_wf,OptWF)
       end
     end
     close(fid_record)
     close(fid_wf)

   elseif Opt["V"] == false && Opt["P"] == false && Opt["Txz"] == true

     #========= write Recordings info =========#
     fid_recinfo = open(path_recinfo, "a+")
     write(fid_recinfo,1803)
     write(fid_recinfo,rn) # Total Receiver Number
     write(fid_recinfo,Medium.Tn) # Total Time Sampling Tn
     write(fid_recinfo,receivers)
     close(fid_recinfo)
     #==========================================#

     #========= write Wavefields info ===========#
     fid_wfinfo = open(path_wfinfo, "a+")
     write(fid_wfinfo,1803)
     write(fid_wfinfo,WFSize.N_BDTpp)
     write(fid_wfinfo,length(Slices))
     write(fid_wfinfo,Slices.*Medium.dt)
     close(fid_wfinfo)
     #==========================================#
     RecLoc_Txz = (PML_Rec[2,:]-1).*WFSize.N_BDTxz[1] + PML_Rec[1,:]
     fid_record = open(path_rec,"a+")
     fid_wf = open(path_wf,"a+")
     for it = 1 : Tn
       Addsource!(WF1, source, it)
       Onestep2D!(WF1, WF2, FD, pml)
       OptRec = WF1.VecBDTxz[RecLoc_Txz]
       write(fid_record,OptRec)
       if it in Slices
         OptWF = WF1.VecBDTxz
         write(fid_wf,OptWF)
       end
     end
     close(fid_record)
     close(fid_wf)

   elseif Opt["V"] == true && Opt["P"] == true && Opt["Txz"] == true

     #========= write Recordings info =========#
     fid_recinfo = open(path_recinfo, "a+")
     write(fid_recinfo,1807)
     write(fid_recinfo,rn) # Total Receiver Number
     write(fid_recinfo,Medium.Tn) # Total Time Sampling Tn
     write(fid_recinfo,receivers)
     close(fid_recinfo)
     #==========================================#

     #========= write Wavefields info ===========#
     fid_wfinfo = open(path_wfinfo, "a+")
     write(fid_wfinfo,1807)
     write(fid_wfinfo,WFSize.N_BDTpp)
     write(fid_wfinfo,length(Slices))
     write(fid_wfinfo,Slices.*Medium.dt)
     close(fid_wfinfo)
     #==========================================#
     RecLoc_Vx = (PML_Rec[2,:]-1).*WFSize.N_BDVx[1] + PML_Rec[1,:]
     RecLoc_Vz = (PML_Rec[2,:]-1).*WFSize.N_BDVz[1] + PML_Rec[1,:]
     RecLoc_P = (PML_Rec[2,:]-1).*WFSize.N_BDTpp[1] + PML_Rec[1,:]
     RecLoc_Txz = (PML_Rec[2,:]-1).*WFSize.N_BDTxz[1] + PML_Rec[1,:]
     fid_record = open(path_rec,"a+")
     fid_wf = open(path_wf,"a+")
     for it = 1 : Tn
       Addsource!(WF1, source, it)
       Onestep2D!(WF1, WF2, FD, pml)
       OptRec = [WF1.VecBDVx[RecLoc_Vx];WF1.VecBDVz[RecLoc_Vz];(WF1.VecBDTxx[RecLoc_P]+WF1.VecBDTzz[RecLoc_P])/2;WF1.VecBDTxz[RecLoc_Txz]]
       write(fid_record,OptRec)
       if it in Slices
         OptWF = [WF1.VecBDVx;WF1.VecBDVz;(WF1.VecBDTxx+WF1.VecBDTzz)/2;WF1.VecBDTxz]
         write(fid_wf,OptWF)
       end
     end
     close(fid_record)
     close(fid_wf)

   elseif Opt["V"] == true && Opt["P"] == true && Opt["Txz"] == false

     #========= write Recordings info =========#
     fid_recinfo = open(path_recinfo, "a+")
     write(fid_recinfo,1804)
     write(fid_recinfo,rn) # Total Receiver Number
     write(fid_recinfo,Medium.Tn) # Total Time Sampling Tn
     write(fid_recinfo,receivers)
     close(fid_recinfo)
     #==========================================#

     #========= write Wavefields info===========#
     fid_wfinfo = open(path_wfinfo, "a+")
     write(fid_wfinfo,1804)
     write(fid_wfinfo,WFSize.N_BDTpp)
     write(fid_wfinfo,length(Slices))
     write(fid_wfinfo,Slices.*Medium.dt)
     close(fid_wfinfo)
     #==========================================#
     RecLoc_Vx = (PML_Rec[2,:]-1).*WFSize.N_BDVx[1] + PML_Rec[1,:]
     RecLoc_Vz = (PML_Rec[2,:]-1).*WFSize.N_BDVz[1] + PML_Rec[1,:]
     RecLoc_P = (PML_Rec[2,:]-1).*WFSize.N_BDTpp[1] + PML_Rec[1,:]
     fid_record = open(path_rec,"a+")
     fid_wf = open(path_wf,"a+")
     for it = 1 : Tn
       Addsource!(WF1, source, it)
       Onestep2D!(WF1, WF2, FD, pml)
       OptRec = [WF1.VecBDVx[RecLoc_Vx];WF1.VecBDVz[RecLoc_Vz];(WF1.VecBDTxx[RecLoc_P]+WF1.VecBDTzz[RecLoc_P])/2]
       write(fid_record,OptRec)
       if it in Slices
         OptWF = [WF1.VecBDVx;WF1.VecBDVz;(WF1.VecBDTxx+WF1.VecBDTzz)/2]
         write(fid_wf,OptWF)
       end
     end
     close(fid_record)
     close(fid_wf)

   elseif Opt["V"] == true && Opt["P"] == false && Opt["Txz"] == true

     #========= write Recordings info =========#
     fid_recinfo = open(path_recinfo, "a+")
     write(fid_recinfo,1805)
     write(fid_recinfo,rn) # Total Receiver Number
     write(fid_recinfo,Medium.Tn) # Total Time Sampling Tn
     write(fid_recinfo,receivers)
     close(fid_recinfo)
     #==========================================#

     #========= write Wavefields info ===========#
     fid_wfinfo = open(path_wfinfo, "a+")
     write(fid_wfinfo,1805)
     write(fid_wfinfo,WFSize.N_BDTpp)
     write(fid_wfinfo,length(Slices))
     write(fid_wfinfo,Slices.*Medium.dt)
     close(fid_wfinfo)
     #==========================================#
     RecLoc_Vx = (PML_Rec[2,:]-1).*WFSize.N_BDVx[1] + PML_Rec[1,:]
     RecLoc_Vz = (PML_Rec[2,:]-1).*WFSize.N_BDVz[1] + PML_Rec[1,:]
     RecLoc_Txz = (PML_Rec[2,:]-1).*WFSize.N_BDTxz[1] + PML_Rec[1,:]
     fid_record = open(path_rec,"a+")
     fid_wf = open(path_wf,"a+")
     for it = 1 : Tn
       Addsource!(WF1, source, it)
       Onestep2D!(WF1, WF2, FD, pml)
       OptRec = [WF1.VecBDVx[RecLoc_Vx];WF1.VecBDVz[RecLoc_Vz];WF1.VecBDTxz[RecLoc_Txz]]
       write(fid_record,OptRec)
       if it in Slices
         OptWF = [WF1.VecBDVx;WF1.VecBDVz;WF1.VecBDTxz]
         write(fid_wf,OptWF)
       end
     end
     close(fid_record)
     close(fid_wf)

   elseif Opt["V"] == false && Opt["P"] == true && Opt["Txz"] == true

     #========= write Recordings info =========#
     fid_recinfo = open(path_recinfo, "a+")
     write(fid_recinfo,1806)
     write(fid_recinfo,rn) # Total Receiver Number
     write(fid_recinfo,Medium.Tn) # Total Time Sampling Tn
     write(fid_recinfo,receivers)
     close(fid_recinfo)
     #==========================================#

     #========= write Wavefields info ===========#
     fid_wfinfo = open(path_wfinfo, "a+")
     write(fid_wfinfo,1806)
     write(fid_wfinfo,WFSize.N_BDTpp)
     write(fid_wfinfo,length(Slices))
     write(fid_wfinfo,Slices.*Medium.dt)
     close(fid_wfinfo)
     #==========================================#
     RecLoc_P = (PML_Rec[2,:]-1).*WFSize.N_BDTpp[1] + PML_Rec[1,:]
     RecLoc_Txz = (PML_Rec[2,:]-1).*WFSize.N_BDTxz[1] + PML_Rec[1,:]
     fid_record = open(path_rec,"a+")
     fid_wf = open(path_wf,"a+")
     for it = 1 : Tn
       Addsource!(WF1, source, it)
       Onestep2D!(WF1, WF2, FD, pml)
       OptRec = [(WF1.VecBDTxx[RecLoc_P]+WF1.VecBDTzz[RecLoc_P])/2;WF1.VecBDTxz[RecLoc_Txz]]
       write(fid_record,OptRec)
       if it in Slices
         OptWF = [(WF1.VecBDTxx+WF1.VecBDTzz)/2;WF1.VecBDTxz]
         write(fid_wf,OptWF)
       end
     end
     close(fid_record)
     close(fid_wf)

   else error("Please choose output component from V, P and Txz")
   end
 end
#================================================================================================#

#================= Simulate Multiple Source with Output Wavefields in mp4(movie)/pdf(image) format ===================#
function Simulation{T<:Real}(PF::T, WF1::WF2D,WF2::WF2D,WFSize::WF2DSize, Medium::Model,
  source::Array{Source}, FD::NSFDMtx, pml::dampCoef,
  sn::Int64, Slices::Union{UnitRange{Int64},StepRange{Int64,Int64},Array{Int64,1}, Int64}, OptId::String,
  OptCpnt::Union{String,Array{String,1},Array{String,1}}; wbox=6, hbox=6, clip=1.0,
  cmap="seismic", aspect="auto", interval=40, vmax = 0.05, vmin=-0.05, ext = 10, iflag = 2,
  infotxtpath = nothing, infobinpath = nothing, Optpath = nothing)

  Writeinfo(WFSize,source,sn,PF,iflag,Medium,infotxtpath,infobinpath,ext)
  Tn = Medium.Tn
  if Optpath == nothing
    Optpath = pwd()
  end
  if OptId == "movie"
    AnimWFs(OptPath, WF1, WF2, WFSize, Medium, source, FD, pml, ext, iflag, sn, Slices, OptCpnt; wbox=wbox, hbox=hbox, cmap=cmap, aspect=aspect, interval=interval, vmax=vmax, vmin=vmax)
  elseif OptId == "pdf"
    PdfWFs(OptPath, WF1, WF2, WFSize, Medium, source, FD, pml, ext, iflag, sn, Slices, OptCpnt; wbox=wbox, hbox=hbox, cmap=cmap, aspect=aspect, interval=interval, vmax=vmax, vmin=vmin)
  end
end
#================================================================================================#

#================= Simulate Single Source with with Output Wavefields in mp4(movie)/pdf(image) format ===================#
function Simulation{T<:Real}(PF::T, WF1::WF2D,WF2::WF2D,WFSize::WF2DSize, Medium::Model,
  source::Source, FD::NSFDMtx, pml::dampCoef,
  Slices::Union{UnitRange{Int64},StepRange{Int64,Int64},Array{Int64,1}, Int64}, OptId::String,
  OptCpnt::Union{String,Array{String,1},Array{String,1}}; wbox=6, hbox=6, clip=1.0,
  cmap="seismic", aspect="auto", interval=40, vmax = 0.05, vmin=-0.05, ext = 10, iflag = 2,
  infotxtpath = nothing, infobinpath = nothing, Optpath = nothing)

  Writeinfo(WFSize,source,PF,iflag,Medium,infotxtpath,infobinpath,ext)
  Tn = Medium.Tn
  if Optpath == nothing
    Optpath = pwd()
  end
  if OptId == "movie"
    AnimWFs(OptPath, WF1, WF2, WFSize, Medium, source, FD, pml, ext, iflag, Slices, OptCpnt; wbox=wbox, hbox=hbox, cmap=cmap, aspect=aspect, interval=interval, vmax=vmax, vmin=vmax)
  elseif OptId == "pdf"
    PdfWFs(OptPath, WF1, WF2, WFSize, Medium, source, FD, pml, ext, iflag, Slices, OptCpnt; wbox=wbox, hbox=hbox, cmap=cmap, aspect=aspect, interval=interval, vmax=vmax, vmin=vmin)
  end
end
#================================================================================================#

#================= Simulate Multiple Source with No Output===================#
function Simulation{T<:Real}(PF::T, WF1::WF2D,WF2::WF2D, Medium::Model,
  source::Array{Source}, FD::NSFDMtx, pml::dampCoef,
  sn::Int64; ext = 10, iflag = 2)

  Tn = Medium.Tn
    for it = 1 : Tn
      Addsource!(WF1, sn, source, it)
      Onestep2D!(WF1, WF2, FD, pml)
    end
end
#================================================================================================#

#================= Simulate Single Source with No Output ===================#
function Simulation{T<:Real}(PF::T, WF1::WF2D,WF2::WF2D, Medium::Model,
  source::Source, FD::NSFDMtx, pml::dampCoef; ext = 10, iflag = 2)

  Tn = Medium.Tn
    for it = 1 : Tn
      Addsource!(WF1, source, it)
      Onestep2D!(WF1, WF2, FD, pml)
    end
end
#================================================================================================#
