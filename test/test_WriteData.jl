# # function test_WriteData()
#
#   include("Project/MODL_UNSPLT_PML/src/Output/WriteData.jl")
#
#   a = spzeros(5,1)
#   b = spzeros(5,1)
#
#   a[1]=1.0
#   b[2]=2.0
#   c =
#   path = "Project/MODL_UNSPLT_PML/test/test_WriteData.bin"
#
#   for i = 1:4
#
#     WriteData(path,"w+",sparse([a[i];b[i]]))
#
#   end
#
# # end
#   fid = open(path,"r")
#   t = read(fid,Float64)
#   close(fid)

path = "Project/MODL_UNSPLT_PML/test/test_WriteData.bin"

fid = open(path,"a+")
write(fid,[1 2; 3 4])
close(fid)

fid = open(path,"r+")
t = read(fid,Int64,4)
close(fid)
