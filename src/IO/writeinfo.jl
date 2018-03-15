
#===  ===#
function writeinfo(
    model::Union{elmod2d,acmod2d},
    sou::Array{MTsource};
    PartextPath="SimuPara.txt",
    ParbinPath="SimuPara.bin")

    sounumber = length(sou)
#=== Simu&OutputInfo.txt ===#
    if PartextPath != nothing
        if PartextPath[end-2:end] != "txt"
            error("Check file extensions. It has to be in 'txt' format. ")
        else
            fid = open(PartextPath,"a+")
            wMT(fid,model)
            wST(fid,sounumber,sou,model.medium.pkf)
            wOPT(fid)
            close(fid)
        end
    end
#=== SimuPara.bin ===#
    # disabled, in garbage folder
    # wMSB(model,sou,sounumber,ParbinPath)
end
#======#





#===  ===#
function writeinfo(
    model::Union{elmod2d,acmod2d},
    sou::Array{MTsource},
    rec::receiver,
    OptCpnt::Array{String};
    PartextPath="SimuPara.txt",
    ParbinPath="SimuPara.bin")

    sounumber = length(sou)
#=== Simu&OutputInfo.txt ===#
    if PartextPath != nothing
        if PartextPath[end-2:end] != "txt"
            error("Check file extensions. It has to be in 'txt' format. ")
        else
            fid = open(PartextPath,"a+")
            wMT(fid,model)
            wST(fid,sounumber,sou,model.medium.pkf)
            wOPT(fid,model,rec,OptCpnt)
            close(fid)
        end
    end
# #=== SimuPara.bin ===#
#     wMSB(model,sou,sounumber,ParbinPath)
end
#======#





#===  ===#
function writeinfo(
    model::Union{elmod2d,acmod2d},
    sou::Array{MTsource},
    Slices::Union{StepRange{Int64,Int64},Array{Int64,1}},
    OptCpnt::Array{String};
    PartextPath="SimuPara.txt",
    ParbinPath="SimuPara.bin")

    sounumber = length(sou)
#=== Simu&OutputInfo.txt ===#
    if PartextPath != nothing
        if PartextPath[end-2:end] != "txt"
            error("Check file extensions. It has to be in 'txt' format. ")
        else
            fid = open(PartextPath,"a+")
            wMT(fid,model)
            wST(fid,sounumber,sou,model.medium.pkf)
            wOPT(fid,model,Slices,OptCpnt)
            close(fid)
        end
    end
# #=== SimuPara.bin ===#
#     wMSB(model,sou,sounumber,ParbinPath)
end
#======#









#===  ===#
function writeinfo(
    model::Union{elmod2d,acmod2d},
    sou::Array{MTsource},
    Slices::Union{StepRange{Int64,Int64},Array{Int64,1}},
    rec::receiver,
    OptCpnt::Array{String};
    PartextPath="SimuPara.txt",
    ParbinPath="SimuPara.bin")

    sounumber = length(sou)
#=== Simu&OutputInfo.txt ===#
    if PartextPath != nothing
        if PartextPath[end-2:end] != "txt"
            error("Check file extensions. It has to be in 'txt' format. ")
        else
            fid = open(PartextPath,"a+")
            wMT(fid,model)
            wST(fid,sounumber,sou,model.medium.pkf)
            wOPT(fid,model,Slices,rec,OptCpnt)
            close(fid)
        end
    end
# #=== SimuPara.bin ===#
#     wMSB(model,sou,sounumber,ParbinPath)
end
#======#








#===  ===#
function writeinfo(
    model::Union{elmod2d},
    sou::Array{MTsource},
    Slices::Union{StepRange{Int64,Int64},Array{Int64,1}};
    PartextPath="SimuPara.txt",
    ParbinPath="SimuPara.bin")

    sounumber = length(sou)
#=== Simu&OutputInfo.txt ===#
    if PartextPath != nothing
        if PartextPath[end-2:end] != "txt"
            error("Check file extensions. It has to be in 'txt' format. ")
        else
            fid = open(PartextPath,"a+")
            wMT(fid,model)
            wST(fid,sounumber,sou,model.medium.pkf)
            wOPT(fid,model,Slices)
            close(fid)
        end
    end
# #=== SimuPara.bin ===#
#     wMSB(model,sou,sounumber,ParbinPath)
end
#======#








#===  ===#
function writeinfo(
    model::Union{acmod2d},
    sou::Array{MTsource},
    Slices::Union{StepRange{Int64,Int64},Array{Int64,1}};
    PartextPath="SimuPara.txt",
    ParbinPath="SimuPara.bin")

    sounumber = length(sou)
#=== Simu&OutputInfo.txt ===#
    if PartextPath != nothing
        if PartextPath[end-2:end] != "txt"
            error("Check file extensions. It has to be in 'txt' format. ")
        else
            fid = open(PartextPath,"a+")
            wMT(fid,model)
            wST(fid,sounumber,sou,model.medium.pkf)
            wOPT(fid,model,Slices)
            close(fid)
        end
    end
# #=== SimuPara.bin ===#
#     wMSB(model,sou,sounumber,ParbinPath)
end
#======#
