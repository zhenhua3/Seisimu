#=== wavefield initialization===#
#=== elastic 2d ===#
function initelwf(nDZ::Int64, nHX::Int64, ext::Int64, iflag::Int64)

        if iflag == 1 # free surface
                return elwf2d(zeros(Float64,(nDZ+ext),(nHX+2*ext)),
                              zeros(Float64,(nDZ+ext),(nHX+2*ext)),
                              zeros(Float64,(nDZ-1+ext),(nHX-1+2*ext)),
                              zeros(Float64,(nDZ+ext),(nHX-1+2*ext)),
                              zeros(Float64,(nDZ-1+ext),(nHX+2*ext)),
                              zeros(Float64,(nDZ+ext),(nHX+2*ext)),
                              zeros(Float64,(nDZ-1+ext),(nHX-1+2*ext)),
                              zeros(Float64,(nDZ+ext),(nHX-1+2*ext)),
                              zeros(Float64,(nDZ-1+ext),(nHX+2*ext)))
        elseif iflag == 2 # unlimited medium
                return elwf2d(zeros(Float64,(nDZ+2*ext),(nHX+2*ext)),
                              zeros(Float64,(nDZ+2*ext),(nHX+2*ext)),
                              zeros(Float64,(nDZ-1+2*ext),(nHX-1+2*ext)),
                              zeros(Float64,(nDZ+2*ext),(nHX-1+2*ext)),
                              zeros(Float64,(nDZ-1+2*ext),(nHX+2*ext)),
                              zeros(Float64,(nDZ+2*ext),(nHX+2*ext)),
                              zeros(Float64,(nDZ-1+2*ext),(nHX-1+2*ext)),
                              zeros(Float64,(nDZ+2*ext),(nHX-1+2*ext)),
                              zeros(Float64,(nDZ-1+2*ext),(nHX+2*ext)))
        else error("Please choose from 1 or 2. 1: free surface; 2: unlimited medium")
        end
end

#=== acoustic 2d ===#
function initacwf(nDZ::Int64, nHX::Int64, ext::Int64, iflag::Int64)

    # This function initialize the wavefield for forward modeling
    # It works for 2D cases

        if iflag == 1 # free surface
                return acwf2d(
                        zeros(Float64,(nDZ+ext),(nHX+2*ext)),
                        zeros(Float64,(nDZ+ext),(nHX+2*ext)),
                        zeros(Float64,(nDZ+ext),(nHX-1+2*ext)),
                        zeros(Float64,(nDZ-1+ext),(nHX+2*ext)),
                        zeros(Float64,(nDZ+ext),(nHX+2*ext)),
                        zeros(Float64,(nDZ+ext),(nHX-1+2*ext)),
                        zeros(Float64,(nDZ-1+ext),(nHX+2*ext)))
        elseif iflag == 2 # unlimited medium
                return acwf2d(
                        zeros(Float64,(nDZ+2*ext),(nHX+2*ext)),
                        zeros(Float64,(nDZ+2*ext),(nHX+2*ext)),
                        zeros(Float64,(nDZ+2*ext),(nHX-1+2*ext)),
                        zeros(Float64,(nDZ-1+2*ext),(nHX+2*ext)),
                        zeros(Float64,(nDZ+2*ext),(nHX+2*ext)),
                        zeros(Float64,(nDZ+2*ext),(nHX-1+2*ext)),
                        zeros(Float64,(nDZ-1+2*ext),(nHX+2*ext)))
        else error("Please choose from 1 or 2. 1: free surface; 2: unlimited medium")
        end
end


#=== elastic 3d ===#
function initelwf(nDZ::Int64, nHX::Int64, nHY::Int64, ext::Int64, iflag::Int64)

        if iflag == 1 # free surface
                return elwf3d(SharedArray{Float32,3}((nDZ-1+ext),(nHX+2*ext),(nHY+2*ext)),
                                 SharedArray{Float32,3}((nDZ+ext),(nHX-1+2*ext),(nHY+2*ext)),
                                 SharedArray{Float32,3}((nDZ+ext),(nHX+2*ext),(nHY-1+2*ext)),
                                 SharedArray{Float32,3}((nDZ+ext),(nHX+2*ext),(nHY+2*ext)),
                                 SharedArray{Float32,3}((nDZ+ext),(nHX+2*ext),(nHY+2*ext)),
                                 SharedArray{Float32,3}((nDZ+ext),(nHX+2*ext),(nHY+2*ext)),
                                 SharedArray{Float32,3}((nDZ-1+ext),(nHX-1+2*ext),(nHY+2*ext)),
                                 SharedArray{Float32,3}((nDZ+ext),(nHX-1+2*ext),(nHY-1+2*ext)),
                                 SharedArray{Float32,3}((nDZ-1+ext),(nHX+2*ext),(nHY-1+2*ext)))
        elseif iflag == 2 # unlimited medium
                return elwf3d(SharedArray{Float32,3}((nDZ-1+2*ext),(nHX+2*ext),(nHY+2*ext)),
                                 SharedArray{Float32,3}((nDZ+2*ext),(nHX-1+2*ext),(nHY+2*ext)),
                                 SharedArray{Float32,3}((nDZ+2*ext),(nHX+2*ext),(nHY-1+2*ext)),
                                 SharedArray{Float32,3}((nDZ+2*ext),(nHX+2*ext),(nHY+2*ext)),
                                 SharedArray{Float32,3}((nDZ+2*ext),(nHX+2*ext),(nHY+2*ext)),
                                 SharedArray{Float32,3}((nDZ+2*ext),(nHX+2*ext),(nHY+2*ext)),
                                 SharedArray{Float32,3}((nDZ-1+2*ext),(nHX-1+2*ext),(nHY+2*ext)),
                                 SharedArray{Float32,3}((nDZ+2*ext),(nHX-1+2*ext),(nHY-1+2*ext)),
                                 SharedArray{Float32,3}((nDZ-1+2*ext),(nHX+2*ext),(nHY-1+2*ext)))
        else error("Please choose from 1 or 2. 1: free surface; 2: unlimited medium")
        end
end







#=== acoustic 3d ===#
function initacwf(nDZ::Int64, nHX::Int64, nHY::Int64, ext::Int64, iflag::Int64)

    # This function initialize the wavefield for forward modeling
    # It works for 2D cases

        if iflag == 1 # free surface
                return acwf3d(SharedArray{Float32,3}((nDZ-1+ext),(nHX+2*ext),(nHY+2*ext)),
                                 SharedArray{Float32,3}((nDZ+ext),(nHX-1+2*ext),(nHY+2*ext)),
                                 SharedArray{Float32,3}((nDZ+ext),(nHX+2*ext),(nHY-1+2*ext)),
                                 SharedArray{Float32,3}((nDZ+ext),(nHX+2*ext),(nHY+2*ext)),
                                 SharedArray{Float32,3}((nDZ+ext),(nHX+2*ext),(nHY+2*ext)),
                                 SharedArray{Float32,3}((nDZ+ext),(nHX+2*ext),(nHY+2*ext)))
        elseif iflag == 2 # unlimited medium

                return acwf3d(SharedArray{Float32,3}((nDZ-1+2*ext),(nHX+2*ext),(nHY+2*ext)),
                                 SharedArray{Float32,3}((nDZ+2*ext),(nHX-1+2*ext),(nHY+2*ext)),
                                 SharedArray{Float32,3}((nDZ+2*ext),(nHX+2*ext),(nHY-1+2*ext)),
                                 SharedArray{Float32,3}((nDZ+2*ext),(nHX+2*ext),(nHY+2*ext)),
                                 SharedArray{Float32,3}((nDZ+2*ext),(nHX+2*ext),(nHY+2*ext)),
                                 SharedArray{Float32,3}((nDZ+2*ext),(nHX+2*ext),(nHY+2*ext)))

        else error("Please choose from 1 or 2. 1: free surface; 2: unlimited medium")
        end
end






#=== size of wavefield initialization ===#
#=== elastic 2d ===#
function initnelwf(nDZ::Int64, nHX::Int64, ext::Int64, iflag::Int64)

  if iflag == 1 # free surface
    return nelwf2d([nDZ+ext,nHX+2*ext],
                          [nDZ-1+ext,nHX-1+2*ext],
                          [nDZ+ext,nHX-1+2*ext],
                          [nDZ-1+ext,nHX+2*ext],
                          [nDZ,nHX],
                          [nDZ-1,nHX-1],
                          [nDZ,nHX-1],
                          [nDZ-1,nHX])

  elseif iflag == 2 # unlimited surface
    return nelwf2d([nDZ+2*ext,nHX+2*ext],
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






#=== acoustic 2d ===#
function initnacwf(nDZ::Int64, nHX::Int64, ext::Int64, iflag::Int64)

  if iflag == 1 # free surface
    return nacwf2d([nDZ+ext,nHX+2*ext],
                          [nDZ+ext,nHX-1+2*ext],
                          [nDZ-1+ext,nHX+2*ext],
                          [nDZ,nHX],
                          [nDZ,nHX-1],
                          [nDZ-1,nHX])

  elseif iflag == 2 # unlimited surface
    return nacwf2d([nDZ+2*ext,nHX+2*ext],
                          [nDZ+2*ext,nHX-1+2*ext],
                          [nDZ-1+2*ext,nHX+2*ext],
                          [nDZ,nHX],
                          [nDZ,nHX-1],
                          [nDZ-1,nHX])
  else error("Please choose from 1 or 2. 1: free surface; 2: unlimited medium")
  end

end





#=== elastic 3d ===#
function initnelwf(nDZ::Int64, nHX::Int64, nHY::Int64, ext::Int64, iflag::Int64)

  if iflag == 1 # free surface
    return nelwf3d([nDZ-1+ext,nHX+2*ext,nHY+2*ext],
                   [nDZ+ext,nHX-1+2*ext,nHY+2*ext],
                   [nDZ+ext,nHX+2*ext,nHY-1+2*ext],
                   [nDZ+ext,nHX+2*ext,nHY+2*ext],
                   [nDZ-1+ext,nHX-1+2*ext,nHY+2*ext],
                   [nDZ+ext,nHX-1+2*ext,nHY-1+2*ext],
                   [nDZ-1+ext,nHX+2*ext,nHY-1+2*ext],
                   [nDZ-1,nHX,nHY],
                   [nDZ,nHX-1,nHY],
                   [nDZ,nHX,nHY-1],
                   [nDZ,nHX,nHY],
                   [nDZ-1,nHX-1,nHY],
                   [nDZ,nHX-1,nHY-1],
                   [nDZ-1,nHX,nHY-1])

  elseif iflag == 2 # unlimited surface
    return nelwf3d([nDZ-1+2*ext,nHX+2*ext,nHY+2*ext],
                   [nDZ+2*ext,nHX-1+2*ext,nHY+2*ext],
                   [nDZ+2*ext,nHX+2*ext,nHY-1+2*ext],
                   [nDZ+2*ext,nHX+2*ext,nHY+2*ext],
                   [nDZ-1+2*ext,nHX-1+2*ext,nHY+2*ext],
                   [nDZ+2*ext,nHX-1+2*ext,nHY-1+2*ext],
                   [nDZ-1+2*ext,nHX+2*ext,nHY-1+2*ext],
                   [nDZ-1,nHX,nHY],
                   [nDZ,nHX-1,nHY],
                   [nDZ,nHX,nHY-1],
                   [nDZ,nHX,nHY],
                   [nDZ-1,nHX-1,nHY],
                   [nDZ,nHX-1,nHY-1],
                   [nDZ-1,nHX,nHY-1])
  else error("Please choose from 1 or 2. 1: free surface; 2: unlimited medium")
  end

end







#=== acoustic 3d ===#
function initnacwf(nDZ::Int64, nHX::Int64, nHY::Int64, ext::Int64, iflag::Int64)

  if iflag == 1 # free surface
    return nacwf3d([nDZ-1+ext,nHX+2*ext,nHY+2*ext],
                   [nDZ+ext,nHX-1+2*ext,nHY+2*ext],
                   [nDZ+ext,nHX+2*ext,nHY-1+2*ext],
                   [nDZ+ext,nHX+2*ext,nHY+2*ext],
                   [nDZ-1,nHX,nHY],
                   [nDZ,nHX-1,nHY],
                   [nDZ,nHX,nHY-1],
                   [nDZ,nHX,nHY])

  elseif iflag == 2 # unlimited surface
    return nacwf3d([nDZ-1+2*ext,nHX+2*ext,nHY+2*ext],
                   [nDZ+2*ext,nHX-1+2*ext,nHY+2*ext],
                   [nDZ+2*ext,nHX+2*ext,nHY-1+2*ext],
                   [nDZ+2*ext,nHX+2*ext,nHY+2*ext],
                   [nDZ-1,nHX,nHY],
                   [nDZ,nHX-1,nHY],
                   [nDZ,nHX,nHY-1],
                   [nDZ,nHX,nHY])
  else error("Please choose from 1 or 2. 1: free surface; 2: unlimited medium")
  end

end
