#=== write Single Source ===#
function wSST(fid::IOStream,Sour::Source,PF::T)
    write(fid,"Source\n") # source
    write(fid,"  Number: 1\n")
    write(fid,"  Location Coordinates [Z,X,Y] (meter):
    [$(Sour.pos[1]),$(Sour.pos[2]),$(Sour.pos[3])]\n")
    write(fid,"  Type: $(Sour.stp)\n")
    write(fid,"  Source Origin Time(s): $(Sour.sot)\n")
    write(fid,"  Discrete Source Origin Time: $(Sour.sotn)\n")
    write(fid,"  Peak Frequency(Hz): $(PF)\n\n")
end
