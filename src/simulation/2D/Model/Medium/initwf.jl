type nspwf2d
  # ===== For unsplit PML boundary simulation ===== #
    txx :: SparseMatrixCSC{Float64, Int64} #Vectoreized with PML Bound
    tzz :: SparseMatrixCSC{Float64, Int64} #Vectoreized with PML Bound
    txz :: SparseMatrixCSC{Float64, Int64} #Vectoreized with PML Bound
    vx :: SparseMatrixCSC{Float64, Int64} #Vectoreized with PML Bound
    vz :: SparseMatrixCSC{Float64, Int64} #Vectoreized with PML Bound
end
type spwf2d
  # ===== For split PML boundary simulation ===== #
    vxx :: SparseMatrixCSC{Float64, Int64}
    vxz :: SparseMatrixCSC{Float64, Int64}
    vzx :: SparseMatrixCSC{Float64, Int64}
    vzz :: SparseMatrixCSC{Float64, Int64}
    txxx :: SparseMatrixCSC{Float64, Int64}
    txxz :: SparseMatrixCSC{Float64, Int64}
    tzzx :: SparseMatrixCSC{Float64, Int64}
    tzzz :: SparseMatrixCSC{Float64, Int64}
    txzx :: SparseMatrixCSC{Float64, Int64}
    txzz :: SparseMatrixCSC{Float64, Int64}
end

type nwf2d
  BDntpp :: Array{Int64,1}
  BDntxz :: Array{Int64,1}
  BDnvx :: Array{Int64,1}
  BDnvz :: Array{Int64,1}
  ntpp :: Array{Int64,1}
  ntxz :: Array{Int64,1}
  nvx :: Array{Int64,1}
  nvz :: Array{Int64,1}
end

function initwf(nDZ::Int64, nHX::Int64, ext::Int64, iflag::Int64, mor::Bool)

    # This function initialize the wavefield for forward modeling
    # It works for 2D cases

        if iflag == 1 # free surface
            if mor == false
                return nspwf2d(spzeros((nDZ+ext),(nHX+2*ext)),
                        spzeros((nDZ+ext),(nHX+2*ext)),
                        spzeros((nDZ-1+ext),(nHX-1+2*ext)),
                        spzeros((nDZ+ext),(nHX-1+2*ext)),
                        spzeros((nDZ-1+ext),(nHX+2*ext)))
            elseif mor == true
                return spwf2d(spzeros((nDZ+ext),(nHX-1+2*ext)),
                        spzeros((nDZ+ext),(nHX-1+2*ext)),
                        spzeros((nDZ-1+ext),(nHX+2*ext)),
                        spzeros((nDZ-1+ext),(nHX+2*ext)),
                        spzeros((nDZ+ext),(nHX+2*ext)),
                        spzeros((nDZ+ext),(nHX+2*ext)),
                        spzeros((nDZ+ext),(nHX+2*ext)),
                        spzeros((nDZ+ext),(nHX+2*ext)),
                        spzeros((nDZ-1+ext),(nHX-1+2*ext)),
                        spzeros((nDZ-1+ext),(nHX-1+2*ext)))
                    end

        elseif iflag == 2 # unlimited medium
            if mor == false
                return nspwf2d(spzeros((nDZ+2*ext),(nHX+2*ext)),
                        spzeros((nDZ+2*ext),(nHX+2*ext)),
                        spzeros((nDZ-1+2*ext),(nHX-1+2*ext)),
                        spzeros((nDZ+2*ext),(nHX-1+2*ext)),
                        spzeros((nDZ-1+2*ext),(nHX+2*ext)))
            elseif mor == true
                return spwf2d(spzeros((nDZ+2*ext),(nHX-1+2*ext)),
                        spzeros((nDZ+2*ext),(nHX-1+2*ext)),
                        spzeros((nDZ-1+2*ext),(nHX+2*ext)),
                        spzeros((nDZ-1+2*ext),(nHX+2*ext)),
                        spzeros((nDZ+2*ext),(nHX+2*ext)),
                        spzeros((nDZ+2*ext),(nHX+2*ext)),
                        spzeros((nDZ+2*ext),(nHX+2*ext)),
                        spzeros((nDZ+2*ext),(nHX+2*ext)),
                        spzeros((nDZ-1+2*ext),(nHX-1+2*ext)),
                        spzeros((nDZ-1+2*ext),(nHX-1+2*ext)))
                    end
        else error("Please choose from 1 or 2. 1: free surface; 2: unlimited medium")
        end
end

function initnwf(nDZ::Int64, nHX::Int64, ext::Int64, iflag::Int64)

  if iflag == 1 # free surface
    return nwf2d([nDZ+ext,nHX+2*ext],
                          [nDZ-1+ext,nHX-1+2*ext],
                          [nDZ+ext,nHX-1+2*ext],
                          [nDZ-1+ext,nHX+2*ext],
                          [nDZ,nHX],
                          [nDZ-1,nHX-1],
                          [nDZ,nHX-1],
                          [nDZ-1,nHX])

  elseif iflag == 2 # unlimited surface
    return nwf2d([nDZ+2*ext,nHX+2*ext],
                          [nDZ-1+2*ext,nHX-1+2*ext],
                          [nDZ+2*ext,nHX-1+2*ext],
                          [nDZ-1+2*ext,nHX+2*ext],
                          [nDZ,nHX],
                          [nDZ-1,nHX-1],
                          [nDZ,nHX-1],
                          [nDZ-1,nHX])
  else error("Please choose from 1 or 2. 1: free surface; 2: unlimited medium")
  end

end
