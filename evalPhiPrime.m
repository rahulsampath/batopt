%Local Derivative
function phiPrime = evalPhiPrime(nodeNum, psi)

if nodeNum == 1
  phiPrime = -0.5;
elseif nodeNum == 2
  phiPrime = 0.5;
else
   error('Invalid NodeNum')
end

