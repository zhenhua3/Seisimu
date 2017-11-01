function TestData()

path = "/home/lzh/Project/STGDFD/test/TestData1.bin"
Nz = 40
Nx = 40

GenerateTestData(Nz,Nx,path)

end


# This function is used to generate test data
function GenerateTestData(Nz :: Int64, Nx :: Int64, Ny :: Int64, path)

    TestData = rand(Nz,Nx,Ny)

    fid = open(path,"w")
    write(fid, TestData)
    close(fid)
end

function GenerateTestData(Nz :: Int64, Nx :: Int64, path)

    TestData = rand(Nz,Nx)

    fid = open(path,"w")
    write(fid, TestData)
    close(fid)
end
