# function Receiver
# Output: discretized receiver locations on the grid
# format: Rec = Array{Int64}(Nr,2),  Rec = [NRz, NRx]
#         Nr : receiver number
#         NRz : Array{Int64}(Nr), discretized x component of receiver
#         NRx : Array{Int64}(Nr), discretized z component of receiver
#
# Input: physical receiver locations, spatial interval dz, dx
# format: location = Array{Float64}(Nr,2), location = [Rz, Rx]
#         Nr : receiver number
#         Rz : Array{Int64}(Nr), z component of receiver
#         Rx : Array{Int64}(Nr), x component of receiver
#         dz : sptial interval
#         dx : spatial interval
#        ext : PML boundary on one side
#      iflag : 1: free surface 2: unlimited boundary


function InitRec{T1,T2,T3<:Real}(rn::Int64,receivers::Array{T1}, dz::T2, dx::T3, ext::Int64, iflag::Int64)

  if rn == 1
   Rec = zeros(Int64,2)
   Rec[1] = Int64(round(receivers[1]./dz))
   Rec[2] = Int64(round(receivers[2]./dx))
   if iflag == 1 # free surface
      PML_Rec = [Rec[1]; Rec[2]+ext]
    elseif iflag == 2 # unlimited space
      PML_Rec = [Rec[1]+ext; Rec[2]+ext]
   end
 elseif rn > 1
   Rec = zeros(Int64,2,rn)
   for i = 1:rn
     Rec[1,i] = Int64(round(receivers[1,i]./dz))
     Rec[2,i] = Int64(round(receivers[2,i]./dx))
   end
   if iflag == 1 # free surface
      PML_Rec = [Rec[1,:]'; Rec[2,:]'+ext]
    elseif iflag == 2 # unlimited space
      PML_Rec = [Rec[1,:]'+ext; Rec[2,:]'+ext]
   end
 end

return Rec, PML_Rec
end
