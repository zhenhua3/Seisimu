type WF2D
  # ===== For unsplit PML boundary simulation ===== #
    VecBDTxx :: SparseMatrixCSC{Float64, Int64} #Vectoreized with PML Bound
    VecBDTzz :: SparseMatrixCSC{Float64, Int64} #Vectoreized with PML Bound
    VecBDTxz :: SparseMatrixCSC{Float64, Int64} #Vectoreized with PML Bound
    VecBDVx :: SparseMatrixCSC{Float64, Int64} #Vectoreized with PML Bound
    VecBDVz :: SparseMatrixCSC{Float64, Int64} #Vectoreized with PML Bound
  # ===== For split PML boundary simulation ===== #
    VecBDVxx :: SparseMatrixCSC{Float64, Int64}
    VecBDVxz :: SparseMatrixCSC{Float64, Int64}
    VecBDVzx :: SparseMatrixCSC{Float64, Int64}
    VecBDVzz :: SparseMatrixCSC{Float64, Int64}
    VecBDTxxx :: SparseMatrixCSC{Float64, Int64}
    VecBDTxxz :: SparseMatrixCSC{Float64, Int64}
    VecBDTzzx :: SparseMatrixCSC{Float64, Int64}
    VecBDTzzz :: SparseMatrixCSC{Float64, Int64}
    VecBDTxzx :: SparseMatrixCSC{Float64, Int64}
    VecBDTxzz :: SparseMatrixCSC{Float64, Int64}
end

type WF2DSize
  N_BDTpp :: Array{Int64,1}
  N_BDTxz :: Array{Int64,1}
  N_BDVx :: Array{Int64,1}
  N_BDVz :: Array{Int64,1}
  N_Tpp :: Array{Int64,1}
  N_Txz :: Array{Int64,1}
  N_Vx :: Array{Int64,1}
  N_Vz :: Array{Int64,1}
end

function InitWF2D(Nz::Int64, Nx::Int64, ext::Int64, iflag::Int64)

    # This function initialize the wavefield for forward modeling
    # It works for 2D cases

        if iflag == 1 # free surface
          InWF2D = WF2D(spzeros((Nz+ext)*(Nx+2*ext),1),
                        spzeros((Nz+ext)*(Nx+2*ext),1),
                        spzeros((Nz-1+ext)*(Nx-1+2*ext),1),
                        spzeros((Nz+ext)*(Nx-1+2*ext),1),
                        spzeros((Nz-1+ext)*(Nx+2*ext),1),
                        spzeros((Nz+ext)*(Nx-1+2*ext),1),
                        spzeros((Nz+ext)*(Nx-1+2*ext),1),
                        spzeros((Nz-1+ext)*(Nx+2*ext),1),
                        spzeros((Nz-1+ext)*(Nx+2*ext),1),
                        spzeros((Nz+ext)*(Nx+2*ext),1),
                        spzeros((Nz+ext)*(Nx+2*ext),1),
                        spzeros((Nz+ext)*(Nx+2*ext),1),
                        spzeros((Nz+ext)*(Nx+2*ext),1),
                        spzeros((Nz-1+ext)*(Nx-1+2*ext),1),
                        spzeros((Nz-1+ext)*(Nx-1+2*ext),1))

        elseif iflag == 2 # unlimited medium
          InWF2D = WF2D(spzeros((Nz+2*ext)*(Nx+2*ext),1),
                        spzeros((Nz+2*ext)*(Nx+2*ext),1),
                        spzeros((Nz-1+2*ext)*(Nx-1+2*ext),1),
                        spzeros((Nz+2*ext)*(Nx-1+2*ext),1),
                        spzeros((Nz-1+2*ext)*(Nx+2*ext),1),
                        spzeros((Nz+2*ext)*(Nx-1+2*ext),1),
                        spzeros((Nz+2*ext)*(Nx-1+2*ext),1),
                        spzeros((Nz-1+2*ext)*(Nx+2*ext),1),
                        spzeros((Nz-1+2*ext)*(Nx+2*ext),1),
                        spzeros((Nz+2*ext)*(Nx+2*ext),1),
                        spzeros((Nz+2*ext)*(Nx+2*ext),1),
                        spzeros((Nz+2*ext)*(Nx+2*ext),1),
                        spzeros((Nz+2*ext)*(Nx+2*ext),1),
                        spzeros((Nz-1+2*ext)*(Nx-1+2*ext),1),
                        spzeros((Nz-1+2*ext)*(Nx-1+2*ext),1))
        else error("Please choose from 1 or 2. 1: free surface; 2: unlimited medium")
        end

        return InWF2D
end

function InitWF2DSize(Nz::Int64, Nx::Int64, ext::Int64, iflag::Int64)

  if iflag == 1 # free surface
    InWF2DSize = WF2DSize([Nz+ext,Nx+2*ext],
                          [Nz-1+ext,Nx-1+2*ext],
                          [Nz+ext,Nx-1+2*ext],
                          [Nz-1+ext,Nx+2*ext],
                          [Nz,Nx],
                          [Nz-1,Nx-1],
                          [Nz,Nx-1],
                          [Nz-1,Nx])
                          
  elseif iflag == 2 # unlimited surface
    InWF2DSize = WF2DSize([Nz+2*ext,Nx+2*ext],
                          [Nz-1+2*ext,Nx-1+2*ext],
                          [Nz+2*ext,Nx-1+2*ext],
                          [Nz-1+2*ext,Nx+2*ext],
                          [Nz,Nx],
                          [Nz-1,Nx-1],
                          [Nz,Nx-1],
                          [Nz-1,Nx])
  else error("Please choose from 1 or 2. 1: free surface; 2: unlimited medium")
  end

  return InWF2DSize
end
