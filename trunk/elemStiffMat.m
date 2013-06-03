%Element Stiffness Matrix
function mat = elemStiffMat()

%2-point quadrature rule
gWts = [1; 1];
gPts = [-(1/sqrt(3)); (1/sqrt(3))];

mat = zeros(2, 2);
for r = 1:2
    for c = 1:2
        for g = 1:2
	    mat(r, c) = mat(r, c) + (gWts(g)*evalPhiPrime(r, gPts(g))*evalPhiPrime(c, gPts(g)));
	end
    end
end

