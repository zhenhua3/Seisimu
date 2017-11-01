function Dx2D(nz::Int64, nx::Int64, Nz::Int64, Nx::Int64, FDCoeff::AbstractArray)

    #::SparseMatrixCSC{Float64,Int64}
    # This function calculate $/partial x vx$
    # The function can automatically determine which one to choose between two differential plans
    # Input : WF WaveField used to calculate its derivative
    #         Nx True model size in X direction. It is useful to determine which differential plan to use
    #         GFDM General Finite-Differential Matrix before truncation

    NFDC = length(FDCoeff)
    # Length of the Finite-Difference Coefficient

    GFDM = spzeros(Nx, Nx+NFDC-1)

    for i = 1 : Nx
        # The reason we dont use NFDC/2 is that NFDC/2 :: Float64 which make it not usable as index
        # Index has to be Int

        GFDM[i , i : i+NFDC-1] = copy( FDCoeff[1 : NFDC] )
        # General Finite-Difference Matrix
        # for further usage we have to trunck the Matrix to a proper size

    end

    if is(nx, Nx) == true # The first way of differential when nx == Nx
        FDM = GFDM[1 : nx-1, div(length(FDCoeff),2) : div(length(FDCoeff),2)+nx-1]
        # Truncate the GFDM to a sproper size

        FDM[1:(div(length(FDCoeff),2)-1),:] = FDM[1:(div(length(FDCoeff),2)-1),:] * 0;
        FDM[(nx-div(length(FDCoeff),2)+1):(nx-1),:] = FDM[(nx-div(length(FDCoeff),2)+1):(nx-1),:] * 0;

        elseif is(nx+1, Nx) == true # The second way of differential when nx != Nx
        FDM = GFDM[1 : nx+1, div(length(FDCoeff),2)+1 : div(length(FDCoeff),2)+nx]

        FDM[1:(div(length(FDCoeff),2)+1-1),:] = FDM[1:(div(length(FDCoeff),2)+1-1),:] * 0;
        FDM[(nx-div(length(FDCoeff),2)+1+1):(nx+1),:] = FDM[(nx-div(length(FDCoeff),2)+1+1):(nx+1),:] * 0;
      else error("Please check your model size")
    end

    I = speye(nz)

    Dx = kron(FDM,I) #kronecker tensor product

    return Dx
end
