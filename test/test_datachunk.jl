using seisimu
Medium_Type = "elastic"
Data_Type = "double"
Nz = 200
Nx = 200
Ny = 208
GPU_mem = 1
CPU_mem = 14
Num_stream = 5

group_full, group_half, stream_full, stream_half = gpuchunk(Medium_Type,Nz,Nx,Ny,Data_Type,GPU_mem,CPU_mem;Num_stream=Num_stream)
