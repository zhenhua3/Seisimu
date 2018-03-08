
#=== non split elastic 2d ===#
function wfpara(medium::elastic2d)

    # Vx ##########################
    AvgRhoVx = zeros(medium.BDnHX,medium.BDnHX-1)
      for i = 1:medium.BDnHX-1
        AvgRhoVx[i:i+1,i] = AvgRhoVx[i:i+1,i] + [0.5,0.5]
      end
    rho_Vx = medium.rho*AvgRhoVx
    rho_Vx = medium.dt./rho_Vx

    #Vz ##########################
    AvgRhoVz = zeros(medium.BDnDZ-1,medium.BDnDZ)
    for i = 1:medium.BDnDZ-1
      AvgRhoVz[i:i,i:i+1] = AvgRhoVz[i:i,i:i+1] .+ [0.5 0.5]
    end
    rho_Vz = AvgRhoVz*medium.rho
    rho_Vz = medium.dt./rho_Vz

    # elastic tpp ######
    lamu = (medium.lambda.+2.*medium.mu).*medium.dt

    #Txz ##########################
    AvgMuTxzX = zeros(medium.BDnHX,medium.BDnHX-1)
    for i = 1:medium.BDnHX-1
      AvgMuTxzX[i:i+1,i] = AvgMuTxzX[i:i+1,i] + [0.25,0.25]
    end
    AvgMuTxzZ = zeros(medium.BDnDZ-1,medium.BDnDZ)
    for i = 1:medium.BDnDZ-1
      AvgMuTxzZ[i:i,i:i+1] = AvgMuTxzZ[i:i,i:i+1] + [1 1]
    end
    mu_Txz = (AvgMuTxzZ*medium.mu*AvgMuTxzX).*medium.dt

    return elwfpara2d(rho_Vz, rho_Vx, lamu, mu_Txz)
end


#=== non split acoustic 2d ===#
function wfpara(medium::acoustic2d)

    # Vx ##########################
    AvgRhoVx = zeros(medium.BDnHX,medium.BDnHX-1)
      for i = 1:medium.BDnHX-1
        AvgRhoVx[i:i+1,i] = AvgRhoVx[i:i+1,i] + [0.5,0.5]
      end
    rho_Vx = medium.rho*AvgRhoVx
    rho_Vx = medium.dt./rho_Vx

    #Vz ##########################
    AvgRhoVz = zeros(medium.BDnDZ-1,medium.BDnDZ)
    for i = 1:medium.BDnDZ-1
      AvgRhoVz[i:i,i:i+1] = AvgRhoVz[i:i,i:i+1] .+ [0.5 0.5]
    end
    rho_Vz = AvgRhoVz*medium.rho
    rho_Vz = medium.dt./rho_Vz

    return acwfpara2d(rho_Vz,rho_Vx)
end




# === non split elastic 3d === #
function wfpara(medium::elastic3d)

    #Vz ##########################
    AvgRhoVz = zeros(medium.BDnDZ-1,medium.BDnDZ)
    for i = 1:medium.BDnDZ-1
        AvgRhoVz[i:i,i:i+1] = AvgRhoVz[i:i,i:i+1] .+ [0.5 0.5]
    end
    rho_Vz = zeros(medium.BDnDZ-1,medium.BDnHX,medium.BDnHY)
    for i in 1:medium.BDnHY
        rho_Vz[:,:,i] = AvgRhoVz*medium.rho[:,:,i]
    end
    rho_Vz = medium.dt./rho_Vz

    # Vx ##########################
    AvgRhoVx = zeros(medium.BDnHX,medium.BDnHX-1)
    for i = 1:medium.BDnHX-1
      AvgRhoVx[i:i+1,i] = AvgRhoVx[i:i+1,i] + [0.5,0.5]
    end
    rho_Vx = zeros(medium.BDnDZ,medium.BDnHX-1,medium.BDnHY)
    for i in 1:medium.BDnHY
        rho_Vx[:,:,i] = medium.rho[:,:,i]*AvgRhoVx
    end
    rho_Vx = medium.dt./rho_Vx

    # Vy #########################
    AvgRhoVy = zeros(medium.BDnHY,medium.BDnHY-1)
    for i = 1:medium.BDnHY-1
      AvgRhoVx[i:i+1,i] = AvgRhoVx[i:i+1,i] + [0.5,0.5]
    end
    rho_Vy = zeros(medium.BDnDZ,medium.BDnHX,medium.BDnHY-1)
    for i in 1:medium.BDnHX
        rho_Vy[:,i,:] = medium.rho[:,i,:]*AvgRhoVy
    end
    rho_Vy = medium.dt./rho_Vy

    #Tpp ##########################
    lamu = (medium.lambda.+2*medium.mu).*medium.dt


     #Txz ##########################
      AvgMuTxzX = zeros(medium.BDnHX,medium.BDnHX-1)
      for i = 1:medium.BDnHX-1
        AvgMuTxzX[i:i+1,i] = AvgMuTxzX[i:i+1,i] + [0.25,0.25]
      end
      AvgMuTxzZ = zeros(medium.BDnDZ-1,medium.BDnDZ)
      for i = 1:medium.BDnDZ-1
        AvgMuTxzZ[i:i,i:i+1] = AvgMuTxzZ[i:i,i:i+1] + [1 1]
      end
      mu_Txz = zeros(medium.BDnDZ-1,medium.BDnHX-1,medium.BDnHY)
      for i in 1:medium.BDnHY
          mu_Txz[:,:,i] = AvgMuTxzZ*medium.mu[:,:,i]*AvgMuTxzX
      end
      mu_Txz = mu_Txz.*medium.dt

        #Tyz ##########################
        AvgMuTyzY = zeros(medium.BDnHY,medium.BDnHY-1)
        for i = 1:medium.BDnHY-1
            AvgMuTyzY[i:i+1,i] = AvgMuTyzY[i:i+1,i] + [0.25,0.25]
        end
        AvgMuTyzZ = zeros(medium.BDnDZ-1,medium.BDnDZ)
        for i = 1:medium.BDnDZ-1
            AvgMuTyzZ[i:i,i:i+1] = AvgMuTyzZ[i:i,i:i+1] + [1 1]
        end
        mu_Tyz = zeros(medium.BDnDZ-1,medium.BDnHX,medium.BDnHY-1)
        for i in 1:medium.BDnHX
            mu_Tyz[:,i,:] = AvgMuTxzZ*medium.mu[:,i,:]*AvgMuTyzY
        end
        mu_Tyz = mu_Tyz.*medium.dt

        #Txy ##########################
        AvgMuTxyY = zeros(medium.BDnHY,medium.BDnHY-1)
        for i = 1:medium.BDnHY-1
            AvgMuTxyY[i:i+1,i] = AvgMuTxyY[i:i+1,i] + [0.25,0.25]
        end
        AvgMuTxyX = zeros(medium.BDnHX-1,medium.BDnHX)
        for i = 1:medium.BDnHX-1
            AvgMuTxyX[i:i,i:i+1] = AvgMuTxyX[i:i,i:i+1] + [1 1]
        end
        mu_Txy = zeros(medium.BDnDZ,medium.BDnHX-1,medium.BDnHY-1)
        for i in 1:medium.BDnDZ
            mu_Txy[i,:,:] = AvgMuTxyX*medium.mu[i,:,:]*AvgMuTxyY
        end
        mu_Txy = mu_Txy.*medium.dt

    return elwfpara3d(rho_Vz,rho_Vx,rho_Vy,lamu,mu_Txz,mu_Txy,mu_Tyz)
end




# === non split acoustic 3d === #
function wfpara(medium::acoustic3d)

    #Vz ##########################
    AvgRhoVz = zeros(medium.BDnDZ-1,medium.BDnDZ)
    for i = 1:medium.BDnDZ-1
        AvgRhoVz[i:i,i:i+1] = AvgRhoVz[i:i,i:i+1] .+ [0.5 0.5]
    end
    rho_Vz = zeros(medium.BDnDZ-1,medium.BDnHX,medium.BDnHY)
    for i in 1:medium.BDnHY
        rho_Vz[:,:,i] = AvgRhoVz*medium.rho[:,:,i]
    end
    rho_Vz = medium.dt./rho_Vz

    # Vx ##########################
    AvgRhoVx = zeros(medium.BDnHX,medium.BDnHX-1)
    for i = 1:medium.BDnHX-1
      AvgRhoVx[i:i+1,i] = AvgRhoVx[i:i+1,i] + [0.5,0.5]
    end
    rho_Vx = zeros(medium.BDnDZ,medium.BDnHX-1,medium.BDnHY)
    for i in 1:medium.BDnHY
        rho_Vx[:,:,i] = medium.rho[:,:,i]*AvgRhoVx
    end
    rho_Vx = medium.dt./rho_Vx

    # Vy #########################
    AvgRhoVy = zeros(medium.BDnHY,medium.BDnHY-1)
    for i = 1:medium.BDnHY-1
      AvgRhoVx[i:i+1,i] = AvgRhoVx[i:i+1,i] + [0.5,0.5]
    end
    rho_Vy = zeros(medium.BDnDZ,medium.BDnHX,medium.BDnHY-1)
    for i in 1:medium.BDnHX
        rho_Vy[:,i,:] = medium.rho[:,i,:]*AvgRhoVy
    end
    rho_Vy = medium.dt./rho_Vy

    return acfd3d(rho_Vz,rho_Vx,rho_Vy)
end
