function extrap2d(model::nspmod2d,
    sou::Array{source};
    rec = nothing,
    RecDataPath = nothing,
    WFDataPath = nothing,
    Slices = nothing,
    OptCpnt = nothing)

    message(rec,RecDataPath,Slices,WFDataPath,OptCpnt)
#=== No output ===#
    if rec == nothing && Slices == nothing

        for it = 1:model.medium.nT
            addsou!(model.wf, sou, it)
            Onestep2D!(model)
            if mod(it,40)==0
                plt.show(plt.imshow(reshape(model.wf.vx,model.nwf.BDnvx[1],model.nwf.BDnvx[2])))
            end
        end


#=== Output Recordings only ===#
  elseif rec != nothing && Slices == nothing
      simurec(model,sou,rec,OptCpnt,RecDataPath)


#=== Output Wavefields only ===#
  elseif rec == nothing && Slices != nothing
      simuwf(model,source,Slices,OptCpnt,WFDataPath)


#=== Output both Recordings and Wavefields ===#
  elseif rec != nothing && Slices != nothing
      simurecwf(model,source,rec,RecDataPath,Slices,WFDataPath,OptCpnt)
  end
end




#=== model order redution ===#
function extrap2d(model::spmod2d,
    source::Array{source};
    WFDataPath = nothing,
    Slices = nothing)
    message(Slices,WFDataPath)
    Extrap2DWF(model,source,rec,RecDataPath)
end



#================= Simulate Multiple Source with Output Wavefields in mp4(movie)/pdf(image) format ===================#
# function Extrap2D(PF::T, WF1::WF2D,WF2::WF2D,WFSize::WF2DSize,
#     Medium::Model,
#   source::Array{Source}, FD::NSFDMtx, pml::dampCoef,
#   sn::Int64,
#   Slices::Union{UnitRange{Int64},StepRange{Int64,Int64},Array{Int64,1}, Int64},
#   OptId::String,
#   OptCpnt::Union{String,Array{String,1},Array{String,1}}; wbox=6, hbox=6, clip=1.0,
#   cmap="seismic", aspect="auto", interval=40, vmax = 0.05,
#   vmin=-0.05, ext = 10, iflag = 2,
#   infotxtpath = nothing, infobinpath = nothing, Optpath = nothing)
#
#   #Writeinfo(WFSize,source,sn,PF,iflag,Medium,infotxtpath,infobinpath,ext)
#   Tn = Medium.Tn
#   if Optpath == nothing
#     Optpath = pwd()
#   end
#   if OptId == "movie"
#     CheckCpnt2D(OptCpnt)
#     sn==1?begin println(sn); AnimWFs(Optpath, WF1, WF2, WFSize, Medium, source, FD, pml, ext, iflag, Slices, OptCpnt; wbox=wbox, hbox=hbox, cmap=cmap, aspect=aspect, interval=interval, vmax=vmax, vmin=vmax) end:
#     begin println(sn);AnimWFs(Optpath, WF1, WF2, WFSize, Medium, source, FD, pml, ext, iflag, sn, Slices, OptCpnt; wbox=wbox, hbox=hbox, cmap=cmap, aspect=aspect, interval=interval, vmax=vmax, vmin=vmax) end
#   elseif OptId == "pdf"
#     CheckCpnt2D(OptCpnt)
#     sn==1?PdfWFs(Optpath, WF1, WF2, WFSize, Medium, source, FD, pml, ext, iflag, Slices, OptCpnt; wbox=wbox, hbox=hbox, cmap=cmap, aspect=aspect, interval=interval, vmax=vmax, vmin=vmin):
#     PdfWFs(Optpath, WF1, WF2, WFSize, Medium, source, FD, pml, ext, iflag, sn, Slices, OptCpnt; wbox=wbox, hbox=hbox, cmap=cmap, aspect=aspect, interval=interval, vmax=vmax, vmin=vmin)
#   end
# end
#================================================================================================#
