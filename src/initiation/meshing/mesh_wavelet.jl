###### wavelet #############

#### Ricker ####
function Ricker(fkp,dt)
    # Input:
    #     fp: peak frequency in Hz
    #     dt: sampling interval in sec
	f0 = fkp/2
	nw=2.2/f0/dt
	nw=2*floor(Int,nw/2)+1
	nc=floor(Int,nw/2)
	w = zeros(nw)
	k= collect(1:nw)
    k=vec(k)
	alpha = (nc-k+1)*f0*dt*pi
	beta=alpha.^2
	w = (1.-beta.*2).*exp.(-beta)
	return w
end
