###### wavelet #############

#### Ricker ####
function Ricker(f0,dt)
    # Input:
    #     f0: central frequency in Hz
    #     dt: sampling interval in sec

	nw=2.2/f0/dt
	nw=2*floor(Int,nw/2)+1
	nc=floor(Int,nw/2)
	w = zeros(nw)
	k= collect(1:nw)
  k=vec(k)
	alpha = (nc-k+1)*f0*dt*pi
	beta=alpha.^2
	w = (1.-beta.*2).*exp(-beta)
	return w
end
