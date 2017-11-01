Mu = zeros(10,10).+1

AvgMuTxzX = zeros(10,9)
for i = 1:9
  AvgMuTxzX[i:i+1,i] = AvgMuTxzX[i:i+1,i] + [0.25,0.25]
end
AvgMuTxzZ = zeros(9,10)
for i = 1:9
  AvgMuTxzZ[i,i:i+1] = AvgMuTxzZ[i,i:i+1] + [1 1]
end
mu_Txz = AvgMuTxzZ*Mu*AvgMuTxzX
