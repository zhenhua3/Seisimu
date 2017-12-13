#=== write Simulation Parameters ===#
function wSPT(fid::IOStream)
    write(fid,"Output files:\n\n")
    write(fid,"SimuPara.bin: ") # SimuPara.bin
    write(fid,"simulation parameters file, read using function: ReadSimuPara  \n\n")
end
