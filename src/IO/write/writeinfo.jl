
#=== Output : MST, MSB ===#
function writeinfo(
    model::Union{nspmod2d,spmod2d},
    sou::Array{source};
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
    wMSB(model,sou,sounumber,ParbinPath)
end
#======#
