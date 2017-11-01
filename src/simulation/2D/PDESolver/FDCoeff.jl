function FDCoeff(order)

    # This function calculate FD coefficient for spatial finite difference
    # Currently we define the coefficient as the one for the 4th order
    # search "Finite difference coefficient" online for different order

    if order == 4
        FDC = [1/24 -9/8 9/8 -1/24] # Finite-Difference Coefficient
        else error("please choose choose 4")
    end
end
