function Cvec = formCvec(M, h, dt, p)

%2-point quadrature rule
gWts = [1; 1];
gPts = [-(1/sqrt(3)); (1/sqrt(3))];

Cvec = zeros(M, 1);

scaling = h/2.0;

nd0 = 1;
x0 = 0;
for e = 1:(M - 1)
    nd1 = nd0 + 1;
    for g = 1:2
        x = localToGlobal(x0, gPts(g), h);
        f = source(x, (p*dt));
        Cvec(nd0) = Cvec(nd0) + (scaling*gWts(g)*f*evalPhi(1, gPts(g)));
        Cvec(nd1) = Cvec(nd1) + (scaling*gWts(g)*f*evalPhi(2, gPts(g)));
    end
    nd0 = nd1;
    x0 = x0 + h;
end


