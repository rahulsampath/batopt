function phi = evalPhi(nodeNum, psi)

if nodeNum == 1
 phi = (1.0 - psi)/2.0;
elseif nodeNum == 2
 phi = (1.0 + psi)/2.0;
else
 error('Invalid nodeNum')
end 


