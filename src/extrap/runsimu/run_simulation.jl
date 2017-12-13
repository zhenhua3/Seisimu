include("in_source.jl")
include("out_recording.jl")
include("out_recwf.jl")
include("out_wavefield.jl")

function runsimu(model::Union{elmod2d,acmod2d},
    sou::Array{source};
    rec = nothing,
    showevery = 100,
    RecDataPath = nothing,
    WFDataPath = nothing,
    Slices = nothing,
    OptCpnt = nothing)

    message(rec,RecDataPath,Slices,WFDataPath,OptCpnt)
#=== No output ===#
    if rec == nothing && Slices == nothing
        for it = 1:model.medium.nT
            addsou!(model.wf, sou, it)
            extrap2d!(model)
            # if mod(it,showevery)==0
            #     plt.show(plt.imshow(model.wf.txx))
            # end
        end


#=== Output Recordings only ===#
  elseif rec != nothing && Slices == nothing
      simurec(model,sou,rec,OptCpnt,RecDataPath)


#=== Output Wavefields only ===#
  elseif rec == nothing && Slices != nothing
      simuwf(model,sou,Slices,OptCpnt,WFDataPath)


#=== Output both Recordings and Wavefields ===#
  elseif rec != nothing && Slices != nothing
      simurecwf(model,sou,rec,RecDataPath,Slices,WFDataPath,OptCpnt)
  end
end
