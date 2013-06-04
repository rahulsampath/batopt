function Bmat = formBmat(elemMmat, M, Ne, h, dt, cp)

Bmat = zeros(M, M);

nd0 = 1;
for i = 1:3
    scaling = 0.5*h*cp(i)/dt;
    for e = 1:(Ne(i))
    nd1 = nd0 + 1;
    Bmat(nd0, nd0) = Bmat(nd0, nd0) + (scaling*elemMmat(1,1));
    Bmat(nd0, nd1) = Bmat(nd0, nd1) + (scaling*elemMmat(1,2));
    Bmat(nd1, nd0) = Bmat(nd1, nd0) + (scaling*elemMmat(2,1));
    Bmat(nd1, nd1) = Bmat(nd1, nd1) + (scaling*elemMmat(2,2));
    nd0 = nd1;
    end
end


