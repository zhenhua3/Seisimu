# function test_Model_func1()
#
#     include("/home/lzh/Project/MODL_UNSPLT_PML/src/ModelSetup/Medium.jl")
#
#     v = Float64[3000, 4000, 1700, 2200]
#     rho = [2.5, 2.5]
#     LN = 4
#     Depth = 2000.0
#     Horiz = 3000.0
#     PF = 120.0
#     ext = 10
#     iflag = 1
#
#     medium = ModInit()
#     medium = model(v, rho, LN, Depth,Horiz, PF, medium, ext, iflag)
#
# end

# function test_Model_func2()

#     include("/home/lzh/Project/MODL_UNSPLT_PML/src/Model.jl")

#     vp = Float64[3000, 4000, 1700, 2200]
#     vs = Float64[1200,2400,4500,2000]
#     rho = [2.5, 2.5, 2.5, 2,5]
#     Depth = [400,300,600,200.0]
#     Horiz = 3000.0
#     PF = 120.0
#     ext = 10
#     iflag = 2

#     medium = ModInit()
#     medium = Model(vp,vs, rho, Depth,Horiz, PF, medium, ext, iflag)

#     VP = reshape(medium.VP, medium.PML_NDep*medium.PML_NHor,1)
#     path = "/home/lzh/Project/MODL_UNSPLT_PML/test/VP.bin"
#     fid = open(path,"w")
#     write(fid,VP)
#     close(fid)

#     VS = reshape(medium.VS, medium.PML_NDep*medium.PML_NHor,1)
#     path = "/home/lzh/Project/MODL_UNSPLT_PML/test/VS.bin"
#     fid = open(path,"w")
#     write(fid,VS)
#     close(fid)

#     Rho = reshape(medium.Rho, medium.PML_NDep*medium.PML_NHor,1)
#     path = "/home/lzh/Project/MODL_UNSPLT_PML/test/Rho.bin"
#     fid = open(path,"w")
#     write(fid,Rho)
#     close(fid)

# end

# function test_Model_func3()

#     include("/home/lzh/Project/MODL_UNSPLT_PML/src/Model.jl")

#     vp_path = "/home/lzh/Project/MODL_UNSPLT_PML/test/VP.bin"
#     vs_path = "/home/lzh/Project/MODL_UNSPLT_PML/test/VS.bin"
#     rho_path = "/home/lzh/Project/MODL_UNSPLT_PML/test/Rho.bin"
#     dz = 1.25
#     dx = 1.25
#     NDep = 1220
#     NHor = 2420
#     PF = 120.0
#     ext = 10
#     iflag = 2

#     medium = ModInit()
#     medium = Model(vp_path ,vs_path, rho_path,dz, dx, NDep, NHor,PF, medium, ext, iflag)

# end

# function test_Model_func4()
#
#     include("/home/lzh/Project/MODL_UNSPLT_PML/src/Model.jl")
#
#     vp_path = "/home/lzh/Project/MODL_UNSPLT_PML/test/VP.bin"
#
#     dz = 1.25
#     dx = 1.25
#     rho = 2.5
#     NDep = 1220
#     NHor = 2420
#     PF = 120.0
#     ext = 10
#     iflag = 1
#
#     medium = ModInit()
#     medium = Model(vp_path, rho,dz, dx, NDep, NHor,PF, medium, ext, iflag)
#
# end

# using PyCall
# @pyimport matplotlib.pyplot as plt
# plt.imshow(medium.VP)
# plt.show()
